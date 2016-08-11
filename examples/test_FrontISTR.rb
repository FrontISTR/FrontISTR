#!/usr/bin/ruby
#
# USAGE
# ruby test_FrontISTR.rb <pathname of fistr> <pathname of data>
#
require 'fileutils'

$fistr = File.join(File.expand_path(ARGV[0]),'fistr')

def create_hecmw_ctrl(mesh,cnt,res=nil,vis=nil)
  base = File.basename(mesh,".*")
  res = base + ".res" if res == nil
  vis = base + ".vis" if vis == nil
  File.open("hecmw_ctrl.dat","w"){|aFile|
    aFile.print <<END_OF_HECMW_CTRL
##
## HEC-MW control file for FrontSTR
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

def exec_test(dirname,mesh,cnt)
  base = File.basename(mesh,".*")
  Dir.chdir(dirname){
    create_hecmw_ctrl(mesh,cnt)
    FileUtils.rm('0.log') if File.exists?('0.log')
    puts $fistr
    system($fistr)
    puts "return value = #{$?.exitstatus}"
    return 1 if $?.exitstatus != 0
    FileUtils.cp('0.log',base+".log")
    res = compare_log(base+".log",base+"_correct.log")
    if res != 0
      puts "#{dirname} #{mesh} #{cnt}"
      return res
    end
  }
  return 0
end

def exec_notest(dirname,mesh,cnt)
  base = File.basename(mesh,".*")
  Dir.chdir(dirname){
    create_hecmw_ctrl(mesh,cnt)
    FileUtils.rm('0.log') if File.exists?('0.log')
    puts $fistr
    system($fistr)
    puts "return value = #{$?.exitstatus}"
    return 1 if $?.exitstatus != 0
    FileUtils.cp('0.log',base+"_correct.log")
  }
  return 0
end

# 4.5412-317 形式を許す
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
      when /Global Summary/
        g = data['Global'] = {}
        while line=aFile.gets
          if line =~/\/\//
            ary = line.chomp.split
            key = ary[0].gsub("//","")
            max = to_float(ary[1])
            min = to_float(ary[2])
            g[key] = [max,min]
          else
            break
          end
        end
      when /@Element/
        e = data['Element'] = {}
        while line=aFile.gets
          if line =~/\/\//
            ary = line.chomp.split
            key = ary[0].gsub("//","")
            max = to_float(ary[1])
            min = to_float(ary[2])
            e[key] = [max,min]
          else
            break
          end
        end
      when /Maximum Temperature/
        g = data['Global'] || data['Global']={}
        g['Temperature'] = []
        ary = line.chomp.split(":")
        g['Temperature'] << ary[1].to_f
        line=aFile.gets
      when /Minimum Temperature/
        g = data['Global']
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
    if( (correct_g[k][0] - v[0]).abs > 1.0e-6 )
      puts "#{itemname} #{k} max value not coincident actual #{v[0]} : correct #{correct_g[k][0]}"
      return 1
    end
    if( (correct_g[k][1] - v[1]).abs > 1.0e-6 )
      puts "#{itemname} #{k} min value not coincident actual #{v[1]} : correct #{correct_g[k][1]}"
      return 1
    end
  }
  return 0
end

def compare_log(actual,correct)
  act_data = read_log(actual)
  correct_data = read_log(correct)
  g = correct_data['Global']
  e = correct_data['Element']
  if act_data['Global']
    res = compare_item(act_data['Global'],g,'Global')
    return res if res != 0
  end
  if act_data['Element']
    res = compare_item(act_data['Element'],e,'@Element')
    return res if res != 0
  end
  return 0
end

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
["static/exA","A361.msh","A300.cnt"],
["static/exA","A362.msh","A300.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2])
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
["static/exB","B731.msh","B731.cnt"],
["static/exB","B741.msh","B741.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2])
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
["static/exC","C731.msh","C700.cnt"],
["static/exC","C741.msh","C700.cnt"],
].each{|param|
	res = exec_test(param[0],param[1],param[2])
	exit res if res != 0
}
else
  meshArray = Dir.glob(File.join(ARGV[1],"*.msh"))
  meshArray.each{|meshfile|
    base = File.basename(meshfile,".msh")
    param = [
      ARGV[1],
      base+".msh",
      base+".cnt"
    ]
    p param
    res = exec_test(param[0],param[1],param[2])
    exit res if res != 0
  }
end
