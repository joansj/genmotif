import sys,os

if len(sys.argv)!=2:
    print 'python '+sys.argv[0]+' <output_binary_filename>'
    sys.exit()

fn_bin=sys.argv[1]
path_src='src/'
path_obj='obj/'
ext_src='.cpp'
compiler='g++'
compiler_includes=[path_src,'../../../eigen/eigen-3-2-4/']
compiler_flags=['-O3','-std=c++0x']
linker_includes=[]
linker_flags=[]

print 'Getting sources...'
fn_src=[]
for root,dirs,files in os.walk(path_src):
    for file in files:
        if file.endswith(ext_src):
            fn_src.append(file[:-4])

print 'Erasing previous stuff...'
os.system('rm -rf '+path_obj)
os.system('mkdir '+path_obj)
os.system('rm -f '+fn_bin)

print 'Compiling...'
for fnbase in fn_src:
    call=compiler
    for flag in compiler_flags: call+=' '+flag
    for inc in compiler_includes: call+=' -I '+inc
    call+=' -c -o '+path_obj+fnbase+'.o'
    call+=' '+path_src+fnbase+ext_src
    print call
    os.system(call)
    if not os.path.exists(path_obj+fnbase+'.o'): sys.exit()

print 'Linking...'
call=compiler
for flag in linker_flags: call+=' '+flag
for inc in linker_includes: call+=' -I '+inc
call+=' -o '+fn_bin
for fnbase in fn_src: call+=' '+path_obj+fnbase+'.o'
print call
os.system(call)
if not os.path.exists(fn_bin): sys.exit()

print 'Done!'
