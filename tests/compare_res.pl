#!/usr/bin/perl

my $path_ref  = shift;
my $path_test = shift;
use Data::Dumper;

my $eps = 1e-8;

my %r=func($path_ref);
#print Dumper {%r};
my %t=func($path_test);
#print Dumper {%t};
my $return=0;

foreach my $key(sort keys(%r)){
  my $res;
  if ($r{$key} < $eps){
    $res = $r{$key} - $t{$key};
  }else{
    $res = ($r{$key} - $t{$key})/$r{$key};
  }
  next if (exists $r{"EIGENVALUE"} and $key ne "EIGENVALUE");
  if (abs($res) < 1e-4){
    print "    $key: \n";
    print "      result:   success\n";
    print "      residual: $res   \n";
  }else{
    print "    $key: \n";
    print "      result:   failure\n";
    print "      residual: $res   \n";
    $return=1;
  }
}
exit($return);



sub func{
  my %ref;
  my $node, $elem, $node_header_cnt, $elem_header_cnt, $global_header_cnt;
  my @global_header_num,  @global_header_name;
  my @node_header_num,  @elem_header_num;
  my @node_header_name, @elem_header_name;
  my %norm;
  $path = shift;
  if (! -e $path) { print "    error: FILE_NOT_FOUND\n"; exit(127); }
  open(DF, "< $path") or die("Error");
  while(my $lin = <DF>){
    if ($lin =~ /\*global/) {
      ($global_header_cnt)=split(/\s+/,<DF>);
      @global_header_num=(); @global_header_name=();
      for (my $i=0;$i<(int(($global_header_cnt-1)/10)+1);$i++){ push(@global_header_num,split(/\s+/,<DF>)) }
      for (my $i=0;$i<$global_header_cnt;$i++) { chomp($line=<DF>); push(@global_header_name,$line);}
      my $ii;
      my $pline; $pline += $_ for @global_header_num; $pline = int(($pline+(5-1))/5);
      for (my $j=0;$j<$pline;$j++) { my $line = <DF>; $l .= $line; }
      my @p = split(/\s+/,$l);
      for (my $k=0;$k<$global_header_num[$j];$k++){
        $norm{$global_header_name[$j]} += $p[$ii]**2;
        $ii++;
      }
    }

    if ($lin =~ /\*data/) {
      ($node,$elem)          =split(/\s+/,<DF>);
      ($node_header_cnt,$elem_header_cnt)=split(/\s+/,<DF>);
      #print "$node, $elem\n";
      #print "$node_header_cnt, $elem_header_cnt\n\n";
      @node_header_num=(); @node_header_name=();
      for (my $i=0;$i<(int(($node_header_cnt-1)/10)+1);$i++){ push(@node_header_num,split(/\s+/,<DF>)) }
      for (my $i=0;$i<$node_header_cnt;$i++) { chomp($line=<DF>); push(@node_header_name,$line);}
      #for (my $i=0;$i<$node_header_cnt;$i++) { print $node_header_num[$i].", ";print $node_header_name[$i]."\n";}
      my $pline; $pline += $_ for @node_header_num; $pline = int(($pline+(5-1))/5);
      #print $pline ."\n\n";
      for (my $i=0;$i<$node;$i++) {
        my $nodei=<DF>;
        my $l;
        for (my $j=0;$j<$pline;$j++) { my $line = <DF>; $l .= $line; }
        my @p = split(/\s+/,$l);
        my $ii;
        for (my $j=0;$j<=$#node_header_num;$j++){
          #print $node_header_name[$j]. " " . $node_header_num[$j]."\n";
          for (my $k=0;$k<$node_header_num[$j];$k++){
            #print $p[$ii] . "\n";
            $norm{$node_header_name[$j]} += $p[$ii]**2;
            $ii++;
          }
        }
      }
      @elem_header_num=(); @elem_header_name=();
      for (my $i=0;$i<(int(($elem_header_cnt-1)/10)+1);$i++){ push(@elem_header_num,split(/\s+/,<DF>)) }
      for (my $i=0;$i<$elem_header_cnt;$i++) { chomp($line=<DF>); push(@elem_header_name,$line);}
      #for (my $i=0;$i<$elem_header_cnt;$i++) { print $elem_header_num[$i].", ";print $elem_header_name[$i]."\n";}
      my $pline; $pline += $_ for @elem_header_num; $pline = int(($pline+(5-1))/5);
      for (my $i=0;$i<$elem;$i++) {
        my $elemi=<DF>;
        my $l;
        for (my $j=0;$j<$pline;$j++) { my $line = <DF>; $l .= $line; }
        my @p = split(/\s+/,$l);
        my $ii;
        for (my $j=0;$j<=$#elem_header_num;$j++){
          #print $elem_header_name[$j]. " " . $elem_header_num[$j]."\n";
          for (my $k=0;$k<$elem_header_num[$j];$k++){
            #print $p[$ii] . "\n";
            $norm{$elem_header_name[$j]} += $p[$ii]**2;
            $ii++;
          }
        }
      }
    }
  }
  close(DF);
  foreach my $key(keys(%norm)){
    $norm{$key} = $norm{$key}**0.5;
  }
  return %norm;
}
