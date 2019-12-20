#!/usr/bin/ruby
#
# USAGE
# ruby test_FrontISTR.rb <pathname of fistr> <pathname of data>
#
require 'fileutils'

$part = File.join(File.expand_path(ARGV[0]),'/hecmw1/tools/hecmw_part1')
$fistr = File.join(File.expand_path(ARGV[0]),'/fistr1/fistr1')
$threshold = 1.0e-4
$np = ARGV[2]
$nt = ARGV[3]
if $nt.nil? then
  $nt = "1"
end


#
# mesh と cnt ファイルを与えて自動的に hecmw_ctrl.dat を生成する
#
def create_hecmw_ctrl(mesh,cnt,res=nil,vis=nil)
  base = File.basename(mesh,".*")
  res = base + ".res" if res == nil
  vis = base + ".vis" if vis == nil
  if $np.to_i > 1 then
  File.open("hecmw_ctrl.dat","w"){|aFile|
    aFile.print <<END_OF_HECMW_CTRL
##
## HEC-MW control file for FrontSTR
## Auto created
## #{Time.now.strftime("%Y-%m-%d %H:%M:%S")}
##
!MESH, NAME=fstrMSH,TYPE=HECMW-DIST
#{mesh}
!CONTROL,NAME=fstrCNT
#{cnt}
!RESULT,NAME=fstrRES,IO=OUT
#{res}
!RESULT,NAME=vis_out,IO=OUT
#{vis}
!MESH, NAME=part_in,TYPE=HECMW-ENTIRE
#{mesh}
!MESH, NAME=part_out,TYPE=HECMW-DIST
#{mesh}
END_OF_HECMW_CTRL
}
  File.open("hecmw_part_ctrl.dat","w"){|aFile|
    aFile.print <<END_OF_HECMW_CTRL
!PARTITION,TYPE=NODE-BASED,METHOD=KMETIS,DOMAIN=#{$np}
END_OF_HECMW_CTRL
}
  else
    File.open("hecmw_ctrl.dat","w"){|aFile|
    aFile.print <<END_OF_HECMW_CTRL
##
## HEC-MW control file for FrontISTR
## Auto created
## #{Time.now.strftime("%Y-%m-%d %H:%M:%S")}
##
!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE
#{mesh}
!CONTROL,NAME=fstrCNT
#{cnt}
!RESULT,NAME=fstrRES,IO=OUT
#{res}
!RESULT,NAME=vis_out,IO=OUT
#{vis}
END_OF_HECMW_CTRL
  }

  end

end

#
# FrontISTR をテスト実行する
# 回答の名前を name で与える
# correctLog を与えたらそれと比較、なければ #{name}_correct.log と比較する
# 回答がなければ今回の結果を回答にする
#
def exec_test(dirname,mesh,cnt,name,correctLog=nil)
  Dir.chdir(dirname){
    create_hecmw_ctrl(mesh,cnt)
    print dirname,"/",mesh,", "
    print dirname,"/",cnt,"\n"
    FileUtils.rm('0.log') if File.exists?('0.log')
    if $np.to_i > 1 then
      puts $part
      system($part)
      execcmd = "mpirun --oversubscribe --allow-run-as-root -np " + $np + " " + $fistr + " -t " + $nt
      puts execcmd
      system(execcmd)
    else
      execcmd = $fistr + " -t " + $nt
      puts execcmd
      system(execcmd)
    end

    puts "return value = #{$?.exitstatus}"
    return 1 if $?.exitstatus != 0
    currentLog = name+".log"
    correctLog = name+"_correct.log" if correctLog==nil
    if File.exists? correctLog
      FileUtils.cp('0.log',currentLog)
      res = compare_log(currentLog,correctLog)
      if res != 0
        puts "#{dirname} #{mesh} #{cnt}"
        return res
      end
    else
      FileUtils.cp('0.log',correctLog)
      return -1
    end
  }
  return 0
end

#
# FrontISTR を既存の hecmw_ctrl.dat でテスト実行する
# 回答の名前を name で与える
# 回答がなければ今回の結果を回答にする
#
def exec_test_original(dirname,name)
  Dir.chdir(dirname){
    FileUtils.rm('0.log') if File.exists?('0.log')
    puts $fistr
    system($fistr)
    puts "return value = #{$?.exitstatus}"
    return 1 if $?.exitstatus != 0
    currentLog = name+".log"
    correctLog = name+"_correct.log"
    if File.exists? correctLog
      FileUtils.cp('0.log',currentLog)
      res = compare_log(currentLog,correctLog)
      if res != 0
        puts "#{dirname} #{name}"
        return res
      end
    else
      FileUtils.cp('0.log',correctLog)
      return -1
    end
  }
  return 0
end

# 4.5412-317 形式を許すための特殊 to_float
def to_float(str)
  if str =~ /(\d)([\-\+])/
    str.gsub!($1+$2,$1+"E"+$2)
  end
  return str.to_f
end

def read_log(filename)
  data = {}
  File.open(filename){|aFile|
    line=aFile.gets
    while line
      case(line)
      when /Global Summary :Max\/Min/
        g = data['Node'] = {}
        while line=aFile.gets
          if line =~/\/\//
            ary = line.chomp.split
            key = ary[0].gsub("//","").gsub("13","31")
            max = to_float(ary[1])
            min = to_float(ary[2])
            g[key] = [max,min]
          else
            break
          end
        end
      when /@Element :Max\/Min/
        e = data['Element'] = {}
        while line=aFile.gets
          if line =~/\/\//
            ary = line.chomp.split
            key = ary[0].gsub("//","").gsub("13","31")
            max = to_float(ary[1])
            min = to_float(ary[2])
            e[key] = [max,min]
          else
            break
          end
        end
      when / Global Summary @Node/
        g = data['Node'] = {}
        while line=aFile.gets
          if line =~/\/\//
            ary = line.chomp.split
            key = ary[0].gsub("//","")
            max = to_float(ary[1])
            min = to_float(ary[3])
            g[key] = [max,min]
          else
            break
          end
        end
      when /Global Summary @Element/
        e = data['Element'] = {}
        while line=aFile.gets
          if line =~/\/\//
            ary = line.chomp.split
            key = ary[0].gsub("//","")
            max = to_float(ary[1])
            min = to_float(ary[3])
            e[key] = [max,min]
          else
            break
          end
        end
      when /Maximum Temperature/
        g = data['Node'] || data['Node']={}
        g['Temperature'] = []
        ary = line.chomp.split(":")
        g['Temperature'] << ary[1].to_f
        line=aFile.gets
      when /Minimum Temperature/
        g = data['Node'] || data['Node']={}
        ary = line.chomp.split(":")
        g['Temperature'] << ary[1].to_f
        line=aFile.gets
      else
        line=aFile.gets
      end
    end
  }
  return data
end

def compare_item(actual_g,correct_g,itemname)
  actual_g.each{|k,v|
    if( correct_g.has_key?(k) && (correct_g[k][0] - v[0]).abs > $threshold )
      puts "#{itemname} #{k} max value not coincident actual #{v[0]} : correct #{correct_g[k][0]}"
      return 1
    end
    if( correct_g.has_key?(k) && (correct_g[k][1] - v[1]).abs > $threshold )
      puts "#{itemname} #{k} min value not coincident actual #{v[1]} : correct #{correct_g[k][1]}"
      return 1
    end
  }
  return 0
end

def compare_log(actual,correct)
  act_data = read_log(actual)
  correct_data = read_log(correct)

  g = correct_data['Node']
  e = correct_data['Element']
  if act_data['Node']
    res = compare_item(act_data['Node'],g,'Node')
    return res if res != 0
  end
  if act_data['Element']
    res = compare_item(act_data['Element'],e,'@Element')
    return res if res != 0
  end
  return 0
end

###### main #########

case(ARGV[1])
when("static/exA")
[
["static/exA","A231.msh","A200.cnt"],
["static/exA","A232.msh","A200.cnt"],
["static/exA","A241.msh","A200.cnt"],
["static/exA","A242.msh","A200.cnt"],
["static/exA","A341.msh","A300.cnt"],
["static/exA","A342.msh","A300.cnt"],
["static/exA","A351.msh","A300.cnt"],
["static/exA","A352.msh","A300.cnt"],
["static/exA","A361.msh","A300.cnt","A361_MUMPS_V4_5.log"],
["static/exA","A362.msh","A300.cnt"],
["static/exA","A731.msh","A700.cnt","A731_MUMPS_V4_5.log"],
["static/exA","A741.msh","A700.cnt","A741_MUMPS_V4_5.log"],
["static/exA","A761.msh","A700_33.cnt","A761_MUMPS_V4_5.log"],
["static/exA","A781.msh","A700_33.cnt","A781_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("static/exB")
[
["static/exB","B231.msh","B231.cnt"],
["static/exB","B232.msh","B232.cnt"],
["static/exB","B241.msh","B241.cnt"],
["static/exB","B242.msh","B242.cnt"],
["static/exB","B341.msh","B341.cnt"],
["static/exB","B342.msh","B342.cnt"],
["static/exB","B351.msh","B351.cnt"],
["static/exB","B352.msh","B352.cnt"],
["static/exB","B361.msh","B361.cnt"],
["static/exB","B362.msh","B362.cnt"],
["static/exB","B731.msh","B731.cnt","B731_MUMPS_V4_5.log"],
["static/exB","B741.msh","B741.cnt","B741_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("static/exC")
[
["static/exC","C231.msh","C200.cnt"],
["static/exC","C232.msh","C200.cnt"],
["static/exC","C241.msh","C200.cnt"],
["static/exC","C242.msh","C200.cnt"],
["static/exC","C341.msh","C300.cnt"],
["static/exC","C342.msh","C300.cnt"],
["static/exC","C351.msh","C300.cnt"],
["static/exC","C352.msh","C300.cnt"],
["static/exC","C361.msh","C300.cnt"],
["static/exC","C362.msh","C300.cnt"],
["static/exC","C731.msh","C700.cnt","C731_MUMPS_V4_5.log"],
["static/exC","C741.msh","C700.cnt","C741_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("static/exD")
[
["static/exD","D231.msh","D200.cnt"],
["static/exD","D232.msh","D200.cnt"],
["static/exD","D241.msh","D200.cnt"],
["static/exD","D242.msh","D200.cnt"],
["static/exD","D341.msh","D300.cnt"],
["static/exD","D342.msh","D300.cnt"],
["static/exD","D351.msh","D300.cnt"],
["static/exD","D352.msh","D300.cnt"],
["static/exD","D361.msh","D300.cnt"],
["static/exD","D362.msh","D300.cnt"],
["static/exD","D731.msh","D700.cnt","D731_MUMPS_V4_5.log"],
["static/exD","D741.msh","D700.cnt","D741_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("static/exE")
[
["static/exE","E231.msh","E200.cnt"],
["static/exE","E232.msh","E200.cnt"],
["static/exE","E241.msh","E200.cnt"],
["static/exE","E242.msh","E200.cnt"],
["static/exE","E341.msh","E300.cnt"],
["static/exE","E342.msh","E300.cnt"],
["static/exE","E351.msh","E300.cnt"],
["static/exE","E352.msh","E300.cnt"],
["static/exE","E361.msh","E300.cnt"],
["static/exE","E362.msh","E300.cnt"],
["static/exE","E731.msh","E700.cnt","E731_MUMPS_V4_5.log"],
["static/exE","E741.msh","E700.cnt","E741_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("static/exF")
[
["static/exF","F231.msh","F200.cnt"],
["static/exF","F232.msh","F200.cnt"],
["static/exF","F241.msh","F200.cnt"],
["static/exF","F242.msh","F200.cnt"],
["static/exF","F341.msh","F300.cnt"],
["static/exF","F342.msh","F300.cnt"],
["static/exF","F351.msh","F300.cnt"],
["static/exF","F352.msh","F300.cnt"],
["static/exF","F361.msh","F300.cnt"],
["static/exF","F362.msh","F300.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("static/exG")
[
["static/exG","A231.msh","G200.cnt","A231_MUMPS_V4_5.log"],
["static/exG","A232.msh","G200.cnt","A232_MUMPS_V4_5.log"],
["static/exG","A241.msh","G200.cnt","A241_MUMPS_V4_5.log"],
["static/exG","A242.msh","G200.cnt","A242_MUMPS_V4_5.log"],
["static/exG","A341.msh","G300.cnt","A341_MUMPS_V4_5.log"],
["static/exG","A342.msh","G300.cnt","A342_MUMPS_V4_5.log"],
["static/exG","A351.msh","G300.cnt","A351_MUMPS_V4_5.log"],
["static/exG","A352.msh","G300.cnt","A352_MUMPS_V4_5.log"],
["static/exG","A361.msh","G300.cnt","A361_MUMPS_V4_5.log"],
["static/exG","A362.msh","G300.cnt","A362_MUMPS_V4_5.log"],
["static/exG","A731.msh","G700.cnt","A731_MUMPS_V4_5.log"],
["static/exG","A741.msh","G700.cnt","A741_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("static/exI")
[
["static/exI","A341.msh","I300.cnt"],
["static/exI","A342.msh","I300.cnt"],
["static/exI","A351.msh","I300.cnt"],
["static/exI","A352.msh","I300.cnt"],
["static/exI","A361.msh","I300.cnt"],
["static/exI","A362.msh","I300.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("static/spring_boundary")
	res = exec_test_original("static/spring_boundary","SB")
	exit res if res != 0
when("static/1elem")
[
  ["static/1elem","arruda.msh","arruda.cnt"],
  ["static/1elem","creep.msh","creep.cnt"],
  ["static/1elem","drucker.msh","drucker.cnt"],
  ["static/1elem","mises.msh","mises.cnt"],
  ["static/1elem","mohr.msh","mohr.cnt"],
  ["static/1elem","mohr.msh","mohrshear.cnt"],
  ["static/1elem","neohooke.msh","neohooke.cnt"],
  ["static/1elem","quad001.msh","quad001.cnt"],
  ["static/1elem","ramberg.msh","ramberg.cnt"],
  ["static/1elem","relax.msh","relax.cnt"],
  ["static/1elem","rivlin.msh","rivlin.cnt"],
  ["static/1elem","swift.msh","swift.cnt"],
  ["static/1elem","viscoe.msh","viscoe.cnt"],
  ["static/1elem","viscoe.msh","viscof.cnt"],
].each{|param|
  # この系列は cnt ファイル名を回答名にしていることに注意
  res = exec_test(param[0],param[1],param[2],File.basename(param[2],".*"))
  exit res if res != 0
}
when("eigen/exJ")
[
["eigen/exJ","A231.msh","J200.cnt"],
["eigen/exJ","A232.msh","J200.cnt"],
["eigen/exJ","A241.msh","J200.cnt"],
["eigen/exJ","A242.msh","J200.cnt"],
["eigen/exJ","A341.msh","J300.cnt"],
["eigen/exJ","A342.msh","J300.cnt"],
["eigen/exJ","A351.msh","J300.cnt"],
["eigen/exJ","A352.msh","J300.cnt"],
["eigen/exJ","A361.msh","J300.cnt"],
["eigen/exJ","A362.msh","J300.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("eigen/exK")
[
["eigen/exK","A231.msh","K200.cnt","A231_MUMPS_V4_5.log"],
["eigen/exK","A232.msh","K200.cnt","A232_MUMPS_V4_5.log"],
["eigen/exK","A241.msh","K200.cnt","A241_MUMPS_V4_5.log"],
["eigen/exK","A242.msh","K200.cnt","A242_MUMPS_V4_5.log"],
["eigen/exK","A341.msh","K300.cnt","A341_MUMPS_V4_5.log"],
["eigen/exK","A342.msh","K300.cnt","A342_MUMPS_V4_5.log"],
["eigen/exK","A351.msh","K300.cnt","A351_MUMPS_V4_5.log"],
["eigen/exK","A352.msh","K300.cnt","A352_MUMPS_V4_5.log"],
["eigen/exK","A361.msh","K300.cnt","A361_MUMPS_V4_5.log"],
["eigen/exK","A362.msh","K300.cnt","A362_MUMPS_V4_5.log"],
["eigen/exK","A731.msh","K700.cnt","A731_MUMPS_V4_5.log"],
["eigen/exK","A741.msh","K700.cnt","A741_MUMPS_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"),param[3])
	exit res if res != 0
}
when("heat/exM")
[
["heat/exM","MA361.msh","A.cnt"],
["heat/exM","MB361.msh","B.cnt"],
["heat/exM","MC361.msh","C.cnt"],
["heat/exM","MD361.msh","D.cnt"],
["heat/exM","ME361.msh","E.cnt"],
["heat/exM","MF361.msh","F.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exN")
[
["heat/exN","N231.msh","N.cnt"],
["heat/exN","N232.msh","N.cnt"],
["heat/exN","N241.msh","N.cnt"],
["heat/exN","N242.msh","N.cnt"],
["heat/exN","N341.msh","N.cnt"],
["heat/exN","N342.msh","N.cnt"],
["heat/exN","N351.msh","N.cnt"],
["heat/exN","N352.msh","N.cnt"],
["heat/exN","N361.msh","N.cnt"],
["heat/exN","N362.msh","N.cnt"],
["heat/exN","N731.msh","N.cnt"],
["heat/exN","N741.msh","N.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exO")
[
["heat/exO","O231.msh","O200.cnt"],
["heat/exO","O232.msh","O200.cnt"],
["heat/exO","O241.msh","O200.cnt"],
["heat/exO","O242.msh","O200.cnt"],
["heat/exO","O341.msh","O300.cnt"],
["heat/exO","O342.msh","O300.cnt"],
["heat/exO","O351.msh","O300.cnt"],
["heat/exO","O352.msh","O300.cnt"],
["heat/exO","O361.msh","O300.cnt"],
["heat/exO","O362.msh","O300.cnt"],
["heat/exO","O731.msh","O700.cnt"],
["heat/exO","O741.msh","O700.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exP")
[
["heat/exP","P231.msh","P230.cnt"],
["heat/exP","P232.msh","P230.cnt"],
["heat/exP","P241.msh","P240.cnt"],
["heat/exP","P242.msh","P240.cnt"],
["heat/exP","P341.msh","P340.cnt"],
["heat/exP","P342.msh","P340.cnt"],
["heat/exP","P351.msh","P350.cnt"],
["heat/exP","P352.msh","P350.cnt"],
["heat/exP","P361.msh","P360.cnt"],
["heat/exP","P362.msh","P360.cnt"],
["heat/exP","P731.msh","P700.cnt"],
["heat/exP","P741.msh","P700.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exQ")
[
["heat/exQ","Q231.msh","Q230.cnt"],
["heat/exQ","Q232.msh","Q230.cnt"],
["heat/exQ","Q241.msh","Q240.cnt"],
["heat/exQ","Q242.msh","Q240.cnt"],
["heat/exQ","Q341.msh","Q340.cnt"],
["heat/exQ","Q342.msh","Q340.cnt"],
["heat/exQ","Q351.msh","Q350.cnt"],
["heat/exQ","Q352.msh","Q350.cnt"],
["heat/exQ","Q361.msh","Q360.cnt"],
["heat/exQ","Q362.msh","Q360.cnt"],
["heat/exQ","Q731.msh","Q700.cnt"],
["heat/exQ","Q741.msh","Q700.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exR")
[
["heat/exR","R231.msh","R230.cnt"],
["heat/exR","R232.msh","R230.cnt"],
["heat/exR","R241.msh","R240.cnt"],
["heat/exR","R242.msh","R240.cnt"],
["heat/exR","R341.msh","R340.cnt"],
["heat/exR","R342.msh","R340.cnt"],
["heat/exR","R351.msh","R350.cnt"],
["heat/exR","R352.msh","R350.cnt"],
["heat/exR","R361.msh","R360.cnt"],
["heat/exR","R362.msh","R360.cnt"],
["heat/exR","R731.msh","R700.cnt"],
["heat/exR","R741.msh","R700.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exS")
[
["heat/exS","S231.msh","S.cnt"],
["heat/exS","S232.msh","S.cnt"],
["heat/exS","S241.msh","S.cnt"],
["heat/exS","S242.msh","S.cnt"],
["heat/exS","S341.msh","S.cnt"],
["heat/exS","S342.msh","S.cnt"],
["heat/exS","S351.msh","S.cnt"],
["heat/exS","S352.msh","S.cnt"],
["heat/exS","S361.msh","S.cnt"],
["heat/exS","S362.msh","S.cnt"],
["heat/exS","S731.msh","S.cnt"],
["heat/exS","S741.msh","S.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exT")
[
["heat/exT","T541.msh","T.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exU")
[
["heat/exU","U231.msh","U231.cnt"],
["heat/exU","U232.msh","U232.cnt"],
["heat/exU","U241.msh","U241.cnt"],
["heat/exU","U242.msh","U242.cnt"],
["heat/exU","U341.msh","U341.cnt"],
["heat/exU","U342.msh","U342.cnt"],
["heat/exU","U351.msh","U351.cnt"],
["heat/exU","U352.msh","U352.cnt"],
["heat/exU","U361.msh","U361.cnt"],
["heat/exU","U362.msh","U362.cnt"],
["heat/exU","U731.msh","U731.cnt"],
["heat/exU","U741.msh","U741.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exU2")
[
["heat/exU2","eh1b.msh","eh1b.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("heat/exV")
[
["heat/exV","V342.msh","V342.cnt"],
["heat/exV","V361.msh","V361.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[1],".*"))
	exit res if res != 0
}
when("dynamic/exW")
[
["dynamic/exW","W342_step.msh","W342_c0_ex_m2_t1.cnt"],
["dynamic/exW","W342_step.msh","W342_c0_im_m2_t1.cnt"],
["dynamic/exW","W361_step.msh","W361_c0_ex_m2_t1.cnt"],
["dynamic/exW","W361_step.msh","W361_c0_im_m2_t1.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[2],".*"))
	exit res if res != 0
}
when("dynamic/exX")
[
["dynamic/exX","W342_step.msh","W342_c0_ex_m2_t1.cnt","W342_c0_ex_m2_t1_MUMPS_V4_5.log"],
["dynamic/exX","W342_step.msh","W342_c0_im_m2_t1.cnt","W342_c0_im_m2_t1_MUMPS_V4_5.log"],
["dynamic/exX","W361_step.msh","W361_c0_ex_m2_t1.cnt","W361_c0_ex_m2_t1_MUMPS_V4_5.log"],
["dynamic/exX","W361_step.msh","W361_c0_im_m2_t1.cnt","W361_c0_im_m2_t1_CG_V4_5.log"],
].each{|param|
	res = exec_test(param[0],param[1],param[2],File.basename(param[2],".*"),param[3])
	exit res if res != 0
}
else
	puts "no test implements"
	exit 1
end
