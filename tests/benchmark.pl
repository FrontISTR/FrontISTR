#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use Digest::SHA qw(sha256_hex);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp qw(tempdir);
use Getopt::Long qw(GetOptions);
use JSON::PP;
use POSIX qw(strftime uname);
use Sys::Hostname qw(hostname);
use Time::HiRes qw(time);

my @modes = qw(serial openmp mpi hybrid);
my %generated_cache_keys = map { $_ => 1 } qw(GIT_HASH);

my ($source_dir, $build_dir, $cmake, $build_jobs, $help);
GetOptions(
  'source-dir=s' => \$source_dir,
  'build-dir=s'  => \$build_dir,
  'cmake=s'      => \$cmake,
  'build-jobs=i' => \$build_jobs,
  'help'         => \$help,
) or die "invalid arguments\n";

if ($help) {
  print "Usage: benchmark.pl --source-dir DIR --build-dir DIR --cmake CMD --build-jobs N\n";
  exit 0;
}
die "missing required arguments\n"
  unless defined $source_dir && defined $build_dir && defined $cmake
      && defined $build_jobs;

$source_dir = abs_path($source_dir);
$build_dir  = abs_path($build_dir);
die "source or build directory not found\n" unless $source_dir && $build_dir;

my ($temporary_root, $repository);
my @worktrees;

sub quoted_command {
  return join ' ', map {
    my $arg = $_;
    $arg =~ s/'/'"'"'/g;
    $arg =~ /[^A-Za-z0-9_.,:+=\/@%-]/ ? "'$arg'" : $arg;
  } @_;
}

sub run_command {
  my ($cwd, @command) = @_;
  print '+ ', quoted_command(@command), "\n";
  my $original = getcwd();
  chdir $cwd or die "cannot enter $cwd: $!\n";
  system { $command[0] } @command;
  my $status = $? == -1 ? 127 : $? >> 8;
  chdir $original or die "cannot return to $original: $!\n";
  return $status;
}

sub capture_command {
  my ($cwd, @command) = @_;
  print '+ ', quoted_command(@command), "\n";
  my $original = getcwd();
  chdir $cwd or die "cannot enter $cwd: $!\n";
  open my $pipe, '-|', @command or die "cannot run $command[0]: $!\n";
  local $/;
  my $output = <$pipe> // '';
  close $pipe;
  my $status = $? >> 8;
  chdir $original or die "cannot return to $original: $!\n";
  die "command failed ($status): " . quoted_command(@command) . "\n" if $status;
  $output =~ s/\s+\z//;
  return $output;
}

sub require_success {
  my ($status, @command) = @_;
  die "command failed ($status): " . quoted_command(@command) . "\n" if $status;
}

sub cleanup {
  return unless $temporary_root;
  if ($repository) {
    for my $worktree (reverse @worktrees) {
      run_command($repository, 'git', 'worktree', 'remove', '--force', $worktree);
    }
  }
  remove_tree($temporary_root) if -d $temporary_root;
}

sub read_cache {
  my ($path) = @_;
  open my $file, '<', $path or die "cannot read CMake cache $path: $!\n";
  my %cache;
  while (my $line = <$file>) {
    chomp $line;
    next if $line eq '' || $line =~ m{^(?:#|//)};
    next unless $line =~ /^([^:=]+):([^=]+)=(.*)$/;
    $cache{$1} = { type => $2, value => $3 };
  }
  close $file;
  return \%cache;
}

sub select_configuration {
  my ($cache) = @_;
  my %selected;
  for my $name (sort keys %$cache) {
    next if $cache->{$name}{type} =~ /^(?:INTERNAL|STATIC)$/;
    next if $generated_cache_keys{$name};
    next if $cache->{$name}{value} =~ /-NOTFOUND\z/;
    $selected{$name} = { %{$cache->{$name}} };
  }
  return \%selected;
}

sub extra_cmake_arguments {
  return () unless defined $ENV{BENCHMARK_CMAKE_ARGS};
  require Text::ParseWords;
  return Text::ParseWords::shellwords($ENV{BENCHMARK_CMAKE_ARGS});
}

sub configure_arguments {
  my ($selected) = @_;
  my @arguments;
  for my $name (sort keys %$selected) {
    my $entry = $selected->{$name};
    push @arguments, $entry->{type} eq 'UNINITIALIZED'
      ? "-D$name=$entry->{value}"
      : "-D$name:$entry->{type}=$entry->{value}";
  }
  return (@arguments, extra_cmake_arguments());
}

sub tool_version {
  my ($command) = @_;
  my $output;
  if (open my $pipe, '-|', $command, '--version') {
    $output = <$pipe> // 'unknown';
    close $pipe;
  } else {
    return 'unavailable';
  }
  $output =~ s/[\r\n].*\z//s;
  return $output;
}

sub environment_info {
  my ($selected) = @_;
  my %compilers;
  for my $language (qw(C CXX Fortran)) {
    my $key = "CMAKE_${language}_COMPILER";
    next unless exists $selected->{$key};
    my $path = $selected->{$key}{value};
    $compilers{$language} = { path => $path, version => tool_version($path) };
  }
  my ($system, $node, $release, undef, $machine) = uname();
  my %configuration = map { $_ => $selected->{$_}{value} } keys %$selected;
  my $comparable = {
    system                => $system,
    release               => $release,
    machine               => $machine,
    cmake                 => tool_version($cmake),
    compilers             => \%compilers,
    configuration         => \%configuration,
    extra_cmake_arguments => [extra_cmake_arguments()],
  };
  my $canonical = JSON::PP->new->canonical->encode($comparable);
  return {
    %$comparable,
    hostname    => hostname(),
    fingerprint => 'sha256:' . sha256_hex($canonical),
  };
}

sub configure_and_build {
  my ($source, $build, $selected, $generator) = @_;
  my @configure = ($cmake, '-S', $source, '-B', $build);
  push @configure, ('-G', $generator) if $generator;
  push @configure, configure_arguments($selected);
  require_success(run_command($source, @configure), @configure);

  my @build = ($cmake, '--build', $build, '--parallel', $build_jobs);
  my $build_type = $selected->{CMAKE_BUILD_TYPE}{value};
  push @build, ('--config', $build_type) if defined $build_type && $build_type ne '';
  require_success(run_command($source, @build), @build);
  return read_cache(File::Spec->catfile($build, 'CMakeCache.txt'));
}

sub measure_existing {
  my ($label, $commit, $requested_ref, $build, $selected, $generator, $mode, $only_case, $sample) = @_;
  my $effective_cache = read_cache(File::Spec->catfile($build, 'CMakeCache.txt'));
  $mode //= build_mode($selected);
  die "invalid BENCHMARK_MODE: $mode\n" unless grep { $_ eq $mode } @modes;
  die "$mode requires WITH_MPI=ON\n"
    if $mode =~ /^(mpi|hybrid)$/ && !cache_enabled($selected, 'WITH_MPI');
  die "$mode requires WITH_OPENMP=ON\n"
    if $mode =~ /^(openmp|hybrid)$/ && !cache_enabled($selected, 'WITH_OPENMP');
  my $case = $ENV{BENCHMARK_CASE} // '';
  my $config = read_json("$source_dir/tests/tutorial-benchmark/cases.json");
  my @cases = map { $_->{name} } @{$config->{cases}};
  die "invalid BENCHMARK_CASE: $case\n" if $case && !grep { $_ eq $case } @cases;
  @cases = ($case) if $case;
  @cases = ($only_case) if defined $only_case;
  my $mpi_build = cache_enabled($selected, 'WITH_MPI') ? 1 : 0;
  my %settings = map {
    $_->{name} => { %$_, %{$_->{$mpi_build ? 'mpi' : 'serial'} // {}} }
  } @{$config->{cases}};
  for my $item (@cases) {
    my $dependency = $settings{$item}{requires};
    if ($dependency && !cache_enabled($selected, $dependency)) {
      die "$item requires $dependency=ON\n";
    }
  }
  my $test_build = File::Spec->catdir($build, "tutorial-benchmark-$mode" . ($case ? "-$case" : ''));
  my ($np, $nt) = ($mode =~ /^(mpi|hybrid)$/ ? 2 : 1,
                    $mode =~ /^(openmp|hybrid)$/ ? 2 : 1);
  my (%tests, $run_status);
  print "\n== $label: run $mode tutorial benchmarks ==\n";
  for my $item (@cases) {
    my $setting = $settings{$item};
    my $solver = $setting->{solver} // '';
    my $parameters = join("\n", @{$setting->{solver_parameters} // []});
    my @command = ($cmake, "-DINPUT=$source_dir/tutorial/$item",
      "-DWORK=$test_build/$item" . ($sample ? "/sample-$sample" : ''), "-DSOLVER=$build/fistr1/fistr1",
      "-DPARTITIONER=$build/hecmw1/tools/hecmw_part1", "-DNP=$np", "-DNT=$nt",
      "-DSOLVER_OVERRIDE=$solver", "-DSOLVER_PARAMETERS=$parameters",
      "-DOUTPUT_TYPE=$config->{output_type}",
      '-DMPIEXEC=' . ($selected->{MPIEXEC_EXECUTABLE}{value} // 'mpiexec'),
      '-P', "$source_dir/tests/tutorial-benchmark/run.cmake");
    my $start = time();
    my $status = run_command($source_dir, @command);
    my $seconds = time() - $start;
    $run_status ||= $status;
    $tests{$item} = { duration_seconds => $seconds, mode => $mode, exit_code => $status };
    my $note = $solver;
    $note .= '; parameters: ' . join(' / ', @{$setting->{solver_parameters}}) if $parameters ne '';
    $tests{$item}{solver_override} = $note if $solver ne '';
    $tests{$item}{output_type} = $config->{output_type};
    if ($status) {
      my $log_directory = "$test_build/$item" . ($sample ? "/sample-$sample" : '');
      my $path = "$log_directory/failure.log";
      my @lines;
      if (open my $log, '<', $path) {
        @lines = <$log>;
        close $log;
      }
      my $stage = @lines ? $lines[0] : "Runner failed (exit $status); see preceding error output\n";
      $tests{$item}{error} = @lines > 60
        ? $stage . "[last 60 log lines]\n" . join('', @lines[-60..-1])
        : @lines ? join('', @lines) : $stage;
      $tests{$item}{log_directory} = $log_directory;
    }
    printf "%s %s %s: %.3fs (exit %d)\n", $status ? 'FAILED' : 'PASSED',
      $label, $item, $seconds, $status;
  }
  my %effective = map {
    $_ => $effective_cache->{$_}{value}
  } grep { exists $effective_cache->{$_} } keys %$selected;
  return ({
    requested_ref => $requested_ref,
    commit        => $commit,
    build         => {
      requested       => { map { $_ => $selected->{$_}{value} } keys %$selected },
      extra_arguments => [extra_cmake_arguments()],
      effective       => \%effective,
      generator       => $generator,
    },
    exit_code       => $run_status // 0,
    benchmark_case  => $case,
    run_mode        => $mode,
    tests           => \%tests,
  }, $run_status // 0);
}

sub measure_pair {
  my ($baseline_args, $current_args, $cached) = @_;
  my $config = read_json("$source_dir/tests/tutorial-benchmark/cases.json");
  my @cases = map { $_->{name} } @{$config->{cases}};
  @cases = ($ENV{BENCHMARK_CASE}) if $ENV{BENCHMARK_CASE};
  my ($baseline, $current) = ($cached, undef);
  for my $name (@cases) {
    my @runs;
    my $count = 1;
    for my $sample (1..3) {
      last if $sample > $count;
      for my $side (0, 1) {
        next if $side == 0 && $cached;
        my $args = $side ? $current_args : $baseline_args;
        my ($run) = measure_existing(@$args, $name, $sample);
        $runs[$side] //= { %$run, tests => {} };
        my $test = $run->{tests}{$name};
        my $aggregate = $runs[$side]{tests}{$name} //= { %$test, samples => [] };
        push @{$aggregate->{samples}}, { %$test };
        if ($test->{exit_code}) {
          @{$aggregate}{qw(exit_code error log_directory)} = @{$test}{qw(exit_code error log_directory)};
        }
        if ($side == 0 && $sample == 1) {
          $count = !$test->{exit_code} && $test->{duration_seconds} <= 15 ? 3 : 1;
        }
      }
      if ($cached && $sample == 1) {
        my $test = $cached->{tests}{$name};
        $count = $test && !$test->{exit_code} &&
          (($test->{samples}[0]{duration_seconds} // $test->{duration_seconds}) <= 15) ? 3 : 1;
      }
    }
    for my $side (0, 1) {
      next if $side == 0 && $cached;
      my $test = $runs[$side]{tests}{$name};
      my @times = sort { $a <=> $b } map { $_->{duration_seconds} } @{$test->{samples}};
      $test->{duration_seconds} = $times[int(@times / 2)];
      $test->{unstable} = @times == 3 && $times[-1] - $times[0] >= 0.5
        && $times[-1] - $times[0] >= 0.1 * $test->{duration_seconds} ? JSON::PP::true : JSON::PP::false;
      my $target = $side ? \$current : \$baseline;
      $$target //= { %{$runs[$side]}, tests => {}, exit_code => 0 };
      $$target->{tests}{$name} = $test;
      $$target->{exit_code} ||= $test->{exit_code};
    }
  }
  return ($baseline, $baseline->{exit_code} // 0, $current, $current->{exit_code} // 0);
}

sub classify {
  my ($before, $after, $thresholds, $repeated) = @_;
  $thresholds = { %$thresholds, %{$thresholds->{repeated}} } if $repeated;
  my $seconds = $after - $before;
  my $percent = $before > 0 ? $seconds / $before * 100 : 0;
  my $level = 'unchanged';
  if (abs($seconds) >= $thresholds->{critical_seconds}
      && abs($percent) >= $thresholds->{critical_percent}) {
    $level = 'critical';
  } elsif (abs($seconds) >= $thresholds->{warning_seconds}
           && abs($percent) >= $thresholds->{warning_percent}) {
    $level = 'warning';
  } elsif (abs($seconds) > 0 && abs($percent) >= $thresholds->{note_percent}) {
    $level = 'note';
  }
  $level = 'improvement' if $seconds < 0 && $level ne 'unchanged';
  return ($level, $seconds, $percent);
}

sub compare_runs {
  my ($baseline, $current, $thresholds) = @_;
  my ($before, $after) = ($baseline->{tests}, $current->{tests});
  my @common = sort grep { exists $after->{$_} } keys %$before;
  my @failures;
  for my $entry (['baseline', $baseline], ['current', $current]) {
    my ($label, $run) = @$entry;
    for my $name (sort keys %{$run->{tests}}) {
      my $test = $run->{tests}{$name};
      next unless $test->{exit_code};
      push @failures, { name => $name, revision => $label, commit => $run->{commit},
        exit_code => $test->{exit_code}, error => $test->{error} // 'Execution failed',
        log_directory => $test->{log_directory} };
    }
  }
  my @comparable = grep { !$before->{$_}{exit_code} && !$after->{$_}{exit_code} } @common;
  my %counts = map { $_ => 0 } qw(critical warning note improvement unchanged);
  my %comparisons;
  for my $name (@comparable) {
    my ($level, $seconds, $percent) = classify(
      $before->{$name}{duration_seconds}, $after->{$name}{duration_seconds}, $thresholds,
      @{$before->{$name}{samples} // []} == 3 && @{$after->{$name}{samples} // []} == 3
    );
    $counts{$level}++;
    $comparisons{$name} = {
      baseline_seconds => $before->{$name}{duration_seconds},
      current_seconds  => $after->{$name}{duration_seconds},
      change_seconds   => $seconds,
      change_percent   => $percent,
      severity         => $level,
      mode             => $after->{$name}{mode},
    };
  }

  my %mode_summary;
  for my $mode (@modes) {
    my @names = grep { ($after->{$_}{mode} // '') eq $mode } @comparable;
    my $old_total = 0;
    my $new_total = 0;
    $old_total += $before->{$_}{duration_seconds} for @names;
    $new_total += $after->{$_}{duration_seconds} for @names;
    my ($level, $seconds, $percent) = classify($old_total, $new_total, $thresholds);
    $mode_summary{$mode} = {
      tests            => scalar @names,
      baseline_seconds => $old_total,
      current_seconds  => $new_total,
      change_seconds   => $seconds,
      change_percent   => $percent,
      severity         => $level,
    };
  }

  my %config_diff;
  for my $name (keys %{$current->{build}{effective}}) {
    my $old = $baseline->{build}{effective}{$name};
    my $new = $current->{build}{effective}{$name};
    $config_diff{$name} = { baseline => $old, current => $new }
      if !defined $old || $old ne $new;
  }
  my @added = sort grep { !exists $before->{$_} } keys %$after;
  my @removed = sort grep { !exists $after->{$_} } keys %$before;
  return {
    summary => {
      common  => scalar @common,
      compared => scalar @comparable,
      failed_runs => scalar @failures,
      added   => scalar @added,
      removed => scalar @removed,
      %counts,
    },
    added_tests               => \@added,
    failures                  => \@failures,
    removed_tests             => \@removed,
    configuration_differences => \%config_diff,
    tests                     => \%comparisons,
    modes                     => \%mode_summary,
  };
}

sub print_report {
  my ($comparison, $build_mode, $run_mode) = @_;
  my $summary = $comparison->{summary};
  print "\n== Benchmark comparison ==\n";
  print "build $build_mode, run ", ($run_mode // 'all'), "\n";
  print "common $summary->{common}, added $summary->{added}, removed $summary->{removed}\n";
  print "compared $summary->{compared}, failed runs $summary->{failed_runs}\n";
  if (@{$comparison->{failures}}) {
    print "\n!!! EXECUTION FAILURES (excluded from timings and totals) !!!\n";
    for my $failure (@{$comparison->{failures}}) {
      print "\nFAILED $failure->{revision} $failure->{name} [$failure->{commit}] (exit $failure->{exit_code})\n";
      print "$failure->{error}\n";
      print "Logs: $failure->{log_directory}\n" if $failure->{log_directory};
    }
  }
  print "No successful pairs available for performance comparison.\n" unless $summary->{compared};
  print "WARNING: effective build configurations differ; see JSON\n"
    if keys %{$comparison->{configuration_differences}};
  for my $level_and_heading (
    ['critical', '!!! CRITICAL REGRESSIONS !!!'],
    ['warning',  'WARNINGS'],
  ) {
    my ($level, $heading) = @$level_and_heading;
    my @names = sort {
      $comparison->{tests}{$b}{change_seconds} <=> $comparison->{tests}{$a}{change_seconds}
    } grep { $comparison->{tests}{$_}{severity} eq $level } keys %{$comparison->{tests}};
    next unless @names;
    print "\n$heading (", scalar @names, ")\n";
    for my $name (@names) {
      my $item = $comparison->{tests}{$name};
      printf "  %s: %.3fs -> %.3fs (%+.1f%%, %+.3fs)\n",
        $name, $item->{baseline_seconds}, $item->{current_seconds},
        $item->{change_percent}, $item->{change_seconds};
    }
  }
  print "\nNOTES: $summary->{note} short-running tests changed by at least 20%\n"
    if $summary->{note};
  print "\nParallel-mode totals\n";
  for my $mode (@modes) {
    my $item = $comparison->{modes}{$mode};
    next unless $item->{tests};
    printf "  %-7s %9.3fs -> %9.3fs (%+.1f%%)\n",
      $mode, $item->{baseline_seconds}, $item->{current_seconds}, $item->{change_percent};
  }
}

sub cache_enabled {
  my ($selected, $name) = @_;
  return ($selected->{$name}{value} // '') =~ /^(?:1|ON|TRUE|YES)$/i;
}

sub build_mode {
  my ($selected) = @_;
  return $ENV{BENCHMARK_BUILD_MODE}
    if defined $ENV{BENCHMARK_BUILD_MODE} && length $ENV{BENCHMARK_BUILD_MODE};
  my $mpi = cache_enabled($selected, 'WITH_MPI');
  my $omp = cache_enabled($selected, 'WITH_OPENMP');
  return $mpi && $omp ? 'hybrid' : $mpi ? 'mpi' : $omp ? 'openmp' : 'serial';
}

sub read_json {
  my ($path) = @_;
  open my $file, '<', $path or die "cannot read $path: $!\n";
  local $/;
  my $data = JSON::PP->new->decode(<$file>);
  close $file;
  return $data;
}

sub write_json {
  my ($path, $data) = @_;
  open my $file, '>', $path or die "cannot write $path: $!\n";
  print {$file} JSON::PP->new->canonical->pretty->encode($data);
  close $file;
}

sub finish_report {
  my ($baseline, $current, $thresholds, $output_dir, $build_mode, $run_mode,
      $baseline_status, $current_status) = @_;
  my $comparison = compare_runs($baseline, $current, $thresholds);
  my $report = {
    schema_version => 1,
    generated_at   => strftime('%Y-%m-%dT%H:%M:%SZ', gmtime()),
    build_mode     => $build_mode,
    run_mode       => $run_mode,
    thresholds     => $thresholds,
    baseline       => $baseline,
    current        => $current,
    comparison     => $comparison,
  };
  make_path($output_dir);
  my $stamp = strftime('%Y%m%dT%H%M%SZ', gmtime());
  my $scope = $build_mode . '-' . ($run_mode // 'all');
  my $filename = sprintf '%s-%s-%s-vs-%s.json',
    $stamp, $scope . '-' . ($ENV{BENCHMARK_CASE} // 'tutorials'), substr($current->{commit}, 0, 12),
    substr($baseline->{commit}, 0, 12);
  my $output = File::Spec->catfile($output_dir, $filename);
  write_json($output, $report);

  print_report($comparison, $build_mode, $run_mode);
  print "\nJSON report: $output\n";
  my $failed = $baseline_status || $current_status;
  my $fail_on_regression = lc($ENV{BENCHMARK_FAIL_ON_REGRESSION} // '0');
  $failed = 1 if $fail_on_regression =~ /^(?:1|on|true|yes)$/
    && $comparison->{summary}{critical};
  print "\n== Benchmark summary: ", ($failed ? 'FAIL' : 'PASS'), " ==\n";
  print "build $build_mode / run ", ($run_mode // 'all'), "\n";
  printf "baseline %s / current %s\n", substr($baseline->{commit}, 0, 12),
    substr($current->{commit}, 0, 12);
  printf "%-24s %12s %12s %9s %7s  %s\n", 'Case', 'Baseline', 'Current', 'Change', 'Runs', 'Result';
  my %names = map { $_ => 1 } (keys %{$baseline->{tests}}, keys %{$current->{tests}});
  for my $name (sort keys %names) {
    my ($old, $new) = ($baseline->{tests}{$name}, $current->{tests}{$name});
    my $item = $comparison->{tests}{$name};
    my @times = map {
      !defined $_ ? '-' : $_->{exit_code} ? 'FAILED' : sprintf('%.3fs', $_->{duration_seconds})
    } ($old, $new);
    my $result = $item ? uc($item->{severity})
      : ($old && $old->{exit_code}) || ($new && $new->{exit_code}) ? 'EXECUTION FAILED'
      : 'NOT COMPARED';
    $result = 'WITHIN THRESHOLD' if $result eq 'UNCHANGED';
    $result .= ' [UNSTABLE]' if ($old && $old->{unstable}) || ($new && $new->{unstable});
    my $counts = join('/', map { !defined $_ ? 0 : scalar(@{$_->{samples} // [1]}) } ($old, $new));
    printf "%-24s %12s %12s %9s %7s  %s\n", $name, @times,
      ($item ? sprintf('%+.1f%%', $item->{change_percent}) : '-'), $counts, $result;
  }
  for my $mode (@modes) {
    my $item = $comparison->{modes}{$mode};
    next unless $item->{tests};
    printf "%-24s %11.3fs %11.3fs %+.1f%%\n", "TOTAL ($mode)",
      $item->{baseline_seconds}, $item->{current_seconds}, $item->{change_percent};
  }
  printf "Compared %d cases; failed runs %d; critical %d; warnings %d; notes %d\n",
    @{$comparison->{summary}}{qw(compared failed_runs critical warning note)};
  print "Totals include only cases successful in both revisions.\n";
  print "Times are medians; Runs = baseline/current. UNSTABLE: range >= 10% and >= 0.5s.\n";
  my %output_types;
  for my $run ($baseline, $current) {
    for my $test (values %{$run->{tests}}) {
      $output_types{$test->{output_type}} = 1 if defined $test->{output_type};
    }
  }
  print "OUTPUT_TYPE override: ", join(', ', sort keys %output_types), " (see JSON for each run)\n"
    if %output_types;
  print "Tutorial source files are unchanged.\n";
  for my $name (sort keys %names) {
    my %overrides;
    for my $entry (['baseline', $baseline], ['current', $current]) {
      my ($label, $run) = @$entry;
      next unless exists $run->{tests}{$name};
      my $override = $run->{tests}{$name}{solver_override};
      push @{$overrides{$override}}, $label if defined $override;
    }
    for my $override (sort keys %overrides) {
      print "SOLVER override (", join(' + ', @{$overrides{$override}}), ") $name: $override\n";
    }
  }
  print "Regression failure policy: $fail_on_regression\n";
  return $failed ? 1 : 0;
}

sub load_cached_run {
  my ($path, $commit, $test_suite_fingerprint, $mode) = @_;
  open my $file, '<', $path or die "cannot read $path: $!\n";
  local $/;
  my $data = JSON::PP->new->decode(<$file>);
  close $file;
  for my $key (qw(current baseline run)) {
    next unless ref $data->{$key} eq 'HASH';
    next if ref $data->{$key}{working_tree} eq 'HASH'
      && $data->{$key}{working_tree}{dirty};
    next unless ($data->{$key}{test_suite}{fingerprint} // '')
      eq $test_suite_fingerprint;
    next unless ($data->{$key}{benchmark_case} // '') eq ($ENV{BENCHMARK_CASE} // '');
    next unless ($data->{$key}{run_mode} // '') eq $mode;
    return $data->{$key} if ($data->{$key}{commit} // '') eq $commit;
  }
  die "$path has no measurement for commit $commit\n";
}

sub prepare_baseline {
  my ($directory, $baseline_commit, $reference, $current_commit, $selected,
      $generator, $environment, $current_test_suite, $build_mode) = @_;
  die "baseline preparation requires a clean working tree\n"
    if $current_test_suite->{working_tree};
  my $root = File::Spec->rel2abs($directory, $source_dir);
  die "baseline output already exists: $root\n" if -e $root;
  make_path($root);

  $temporary_root = tempdir('frontistr-benchmark-XXXXXX', TMPDIR => 1, CLEANUP => 0);
  my $worktree = File::Spec->catdir($temporary_root, 'baseline-source');
  my @command = ('git', 'worktree', 'add', '--detach', $worktree, $baseline_commit);
  require_success(run_command($repository, @command), @command);
  push @worktrees, $worktree;

  my $source = File::Spec->catdir($root, 'source');
  my $build = File::Spec->catdir($root, 'build');
  @command = ($cmake, '-E', 'copy_directory', $worktree, $source);
  require_success(run_command($source_dir, @command), @command);
  configure_and_build($source, $build, $selected, $generator);
  write_json(File::Spec->catfile($root, 'metadata.json'), {
    schema_version      => 2,
    build_mode          => $build_mode,
    current_commit      => $current_commit,
    baseline_commit     => $baseline_commit,
    requested_ref       => $reference,
    environment         => $environment,
    current_test_suite  => $current_test_suite,
    baseline_test_suite => $current_test_suite,
  });
  print "\nBaseline build: $root\n";
  return 0;
}

sub compare_prebuilt {
  my ($directory, $current_commit, $selected, $generator, $environment,
      $test_suite, $thresholds, $output_dir, $build_mode, $run_mode) = @_;
  die "BENCHMARK_MODE is required with BENCHMARK_BASELINE_DIR\n"
    unless defined $run_mode && grep { $_ eq $run_mode } @modes;
  my $root = abs_path(File::Spec->rel2abs($directory, $source_dir))
    or die "baseline build not found: $directory\n";
  my $metadata = read_json(File::Spec->catfile($root, 'metadata.json'));
  die "unsupported baseline build metadata\n"
    unless ($metadata->{schema_version} // 0) == 2;
  die "baseline build belongs to a different current commit\n"
    unless ($metadata->{current_commit} // '') eq $current_commit;
  die "baseline build uses a different test suite\n"
    unless ($metadata->{current_test_suite}{fingerprint} // '')
      eq $test_suite->{fingerprint};
  die "baseline build mode does not match $build_mode\n"
    unless ($metadata->{build_mode} // '') eq $build_mode;

  my ($baseline, $baseline_status, $current, $current_status) = measure_pair(
    ['baseline', $metadata->{baseline_commit}, $metadata->{requested_ref},
    File::Spec->catdir($root, 'build'), $selected, $generator, $run_mode],
    ['current', $current_commit, 'HEAD', $build_dir, $selected, $generator, $run_mode]
  );
  $baseline->{environment} = $environment;
  $baseline->{test_suite} = $test_suite;

  $current->{environment} = $environment;
  $current->{test_suite} = $test_suite;
  $current->{working_tree} = {
    dirty       => JSON::PP::false,
    fingerprint => undef,
  };
  return finish_report(
    $baseline, $current, $thresholds, $output_dir, $build_mode, $run_mode,
    $baseline_status, $current_status
  );
}

sub benchmark {
  my $cache = read_cache(File::Spec->catfile($build_dir, 'CMakeCache.txt'));
  my $selected = select_configuration($cache);
  if ($^O eq 'darwin' && exists $selected->{CMAKE_OSX_SYSROOT}
      && -d $selected->{CMAKE_OSX_SYSROOT}{value}) {
    $ENV{SDKROOT} = $selected->{CMAKE_OSX_SYSROOT}{value};
  }
  my $generator = $cache->{CMAKE_GENERATOR}{value};
  my $reference = $ENV{BENCHMARK_REF};
  my $output_dir = $ENV{BENCHMARK_OUTPUT_DIR}
    // File::Spec->catdir($build_dir, 'benchmark-results');
  my $thresholds = {
    note_percent     => 0 + ($ENV{BENCHMARK_NOTE_PERCENT} // 20),
    warning_percent  => 0 + ($ENV{BENCHMARK_WARNING_PERCENT} // 15),
    warning_seconds  => 0 + ($ENV{BENCHMARK_WARNING_SECONDS} // 1),
    critical_percent => 0 + ($ENV{BENCHMARK_CRITICAL_PERCENT} // 10),
    critical_seconds => 0 + ($ENV{BENCHMARK_CRITICAL_SECONDS} // 5),
    repeated => { warning_percent => 5, warning_seconds => 0.5,
      critical_percent => 10, critical_seconds => 1 },
  };

  $repository = capture_command($source_dir, 'git', 'rev-parse', '--show-toplevel');
  my $current_commit = capture_command(
    $repository, 'git', 'rev-parse', '--verify', 'HEAD^{commit}'
  );
  my ($working_tree_dirty, $working_tree_fingerprint) = (0, undef);
  unless (defined $ENV{BENCHMARK_BASELINE_DIR}) {
    my $working_tree_status = capture_command(
      $repository, 'git', 'status', '--short', '--untracked-files=all'
    );
    $working_tree_dirty = $working_tree_status ne '';
    my $untracked_files = capture_command(
      $repository, 'git', 'ls-files', '--others', '--exclude-standard', '-z'
    );
    my @untracked_hashes = map {
      $_ . "\0" . capture_command($repository, 'git', 'hash-object', '--', $_)
    } grep { $_ ne '' } split /\0/, $untracked_files;
    $working_tree_fingerprint = 'sha256:' . sha256_hex(
      $working_tree_status,
      capture_command($repository, 'git', 'diff', '--binary', 'HEAD', '--'),
      @untracked_hashes
    ) if $working_tree_dirty;
  }
  my $test_suite = {
    commit      => $current_commit,
    working_tree => $working_tree_dirty ? JSON::PP::true : JSON::PP::false,
    fingerprint => 'tutorials:' . ($working_tree_fingerprint // "git:$current_commit"),
  };

  my $environment = environment_info($selected);
  my $build_mode = build_mode($selected);
  my $run_mode = $ENV{BENCHMARK_MODE} // $build_mode;
  if (defined $ENV{BENCHMARK_BASELINE_DIR}) {
    return compare_prebuilt(
      $ENV{BENCHMARK_BASELINE_DIR}, $current_commit, $selected, $generator,
      $environment, $test_suite, $thresholds, $output_dir, $build_mode, $run_mode
    );
  }

  $reference = capture_command(
    $repository, 'git', 'describe', '--tags', '--abbrev=0', 'HEAD'
  ) unless defined $reference && length $reference;
  my $baseline_commit = capture_command(
    $repository, 'git', 'rev-parse', '--verify', $reference . '^{commit}'
  );
  my $baseline_test_suite = $test_suite;
  die "current and baseline resolve to the same commit\n"
    if $current_commit eq $baseline_commit;
  if (defined $ENV{BENCHMARK_PREPARE_DIR}) {
    return prepare_baseline(
      $ENV{BENCHMARK_PREPARE_DIR}, $baseline_commit, $reference,
      $current_commit, $selected, $generator, $environment, $test_suite,
      $build_mode
    );
  }

  $temporary_root = tempdir('frontistr-benchmark-XXXXXX', TMPDIR => 1, CLEANUP => 0);
  my ($baseline, $baseline_status, $baseline_args);
  if (defined $ENV{BENCHMARK_BASELINE_JSON}) {
    my $path = abs_path(File::Spec->rel2abs(
      $ENV{BENCHMARK_BASELINE_JSON}, $source_dir
    ))
      or die "baseline JSON not found: $ENV{BENCHMARK_BASELINE_JSON}\n";
    $baseline = load_cached_run(
      $path, $baseline_commit, $baseline_test_suite->{fingerprint}, $run_mode
    );
    die "saved baseline was measured in a different environment\n"
      unless ($baseline->{environment}{fingerprint} // '') eq $environment->{fingerprint};
    print "Using saved baseline from $path\n";
    $baseline_status = $baseline->{exit_code} // 0;
  } else {
    my $source = File::Spec->catdir($temporary_root, 'baseline-source');
    my $build  = File::Spec->catdir($temporary_root, 'baseline-build');
    my @command = ('git', 'worktree', 'add', '--detach', $source, $baseline_commit);
    require_success(run_command($repository, @command), @command);
    push @worktrees, $source;
    configure_and_build($source, $build, $selected, $generator);
    $baseline_args = ['baseline', $baseline_commit, $reference, $build, $selected, $generator, $run_mode];
  }

  my $current_build  = File::Spec->catdir($temporary_root, 'current-build');
  my $current_source = $source_dir;
  my $current_ref = 'working-tree';
  unless ($working_tree_dirty) {
    $current_source = File::Spec->catdir($temporary_root, 'current-source');
    my @command = ('git', 'worktree', 'add', '--detach', $current_source, $current_commit);
    require_success(run_command($repository, @command), @command);
    push @worktrees, $current_source;
    $current_ref = 'HEAD';
  }
  configure_and_build($current_source, $current_build, $selected, $generator);
  my ($current, $current_status);
  ($baseline, $baseline_status, $current, $current_status) = measure_pair(
    $baseline_args,
    ['current', $current_commit, $current_ref, $current_build, $selected, $generator, $run_mode], $baseline
  );
  $baseline->{environment} = $environment;
  $baseline->{test_suite} = $baseline_test_suite;
  $current->{environment} = $environment;
  $current->{test_suite} = { %$test_suite };
  $current->{working_tree} = {
    dirty       => $working_tree_dirty ? JSON::PP::true : JSON::PP::false,
    fingerprint => $working_tree_fingerprint,
  };

  return finish_report(
    $baseline, $current, $thresholds, $output_dir, $build_mode, $run_mode,
    $baseline_status, $current_status
  );
}

my ($status, $error);
eval {
  $status = benchmark();
  1;
} or $error = $@ || "unknown error\n";
cleanup();
if ($error) {
  print STDERR "benchmark: error: $error";
  exit 2;
}
exit $status;
