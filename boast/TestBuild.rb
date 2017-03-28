
require './LibconvBase.rb'
#load './LibconvBase.rb'

dwt=generate_dwt()
iwt=generate_iwt()
mf=generate_mf()
imf=generate_imf()
s1r=generate_s1tor()
rs1=generate_rtos1()

files=[dwt.dump_to_file(),iwt.dump_to_file(),mf.dump_to_file(),imf.dump_to_file(),rs1.dump_to_file(),s1r.dump_to_file()]

raise if files.length() != files.to_set().length()

makelines="""
lib_LIBRARIES = libconv.a

#temporary compiling line for gfortran
AM_FCFLAGS = -I. -O

libconv_a_SOURCES = """

files.each{|f| makelines+=' '+f.to_str}

File::open("src/Makefile.am","w") {|f|
  f.puts makelines
    }







