#!/usr/bin/python3.8

#SBATCH -p low
#SBATCH -N 2
#SBATCH -n 10
#SBATCH -t 3:00:0
#SBATCH -J rizk_pym
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


import os, sys, re, shutil, random, time, pathlib
from pathlib import Path
import numpy as np
import copy
from optparse import OptionParser

starting_directory=os.popen('pwd').read()
starting_directory=starting_directory.rstrip('\n')
print("the starting directory is: \n"+starting_directory)

parser = OptionParser()
parser.add_option("-l", "--ligand", dest="lig", help="ligand name, 3 letters")
parser.add_option("-p", "--pdb", dest="pdb", help="pdb file, path")
parser.add_option("-o", "--origin", dest="origin", default=" ", help="select (pdb name) where to shift the model, leaving this empty will not shift the model")
parser.add_option("-d", "--directory", dest="directory", default="best_solution_1_group", help="directory where to look for merged data")
parser.add_option("-e", "--depth", dest="depth", default="2",help="directory max depth to search in")
parser.add_option("-v", "--setview", dest="view", default=" ", help="add view log file, leaving this empty will take default view of the ligand")
parser.add_option("-c", "--setcolor", dest="color", default=" ",help="ligand color")
parser.add_option("-m", "--mapcolor", dest="map_color", default="tv_blue",help="density map color")
parser.add_option("-t", "--mesh_thickness", dest="thick", default="0.5")
parser.add_option("-f", "--frames", dest="frames", default="60")
parser.add_option("-r", "--refinement_tool", dest="refine", default="refmac",help="refmac or phenix, by default refmac")
parser.add_option("", "--fft-map", dest="fft_map", default="")
parser.add_option("", "--mtz", dest="user_mtz", default="")
parser.add_option("-2", "--rmsd_2fofc", dest="rmsd2fofc", default="1")
parser.add_option("-1", "--rmsd_fofc", dest="rmsdfofc", default="3")
parser.add_option("", "--photo", action="store_true", dest="photo", default=False)
parser.add_option("", "--video", action="store_true", dest="video", default=False)
parser.add_option("", "--gif", action="store_true", dest="gif", default=False)
parser.add_option("", "--full_map", action="store_true", dest="not_omit", default=False)
parser.add_option("", "--hkl", dest="hkl", default="merged.hkl")
parser.add_option("", "--mypath", dest="mypath",default=".", help="path where to look for directories")
parser.add_option("", "--overwrite", action="store_true", dest="overwrite", default=False)

(options, args) = parser.parse_args()


if options.lig is None:
    lig = input ("Enter ligand (or amino acid) name: ")
else:
    lig = options.lig

LIG = lig.upper()

if options.pdb is None:
    pdb = input ("Enter pdb name: ")
else:
    pdb = options.pdb
pdbfullpath=os.path.abspath(pdb)
pdbname=Path(pdb).stem

#                                                          SETTING VIEW
if options.view != " " :
    print('VIEW set by user')
    with open(str(options.view), 'r') as viewlog:
        viewdata = viewlog.read().replace('\n', '')
else:
    viewdata=" "

#                                                          FINDING PATHS
#print(str(options.directory))
paths = []
if (options.directory != ".") and (options.directory != "./"):
    mypath=options.mypath
    orig_path = mypath.count(os.path.sep)
    print("will look for "+str(options.hkl)+" in "+str(options.directory)+" folders at:")
    for mypath, dirs, files, in os.walk(mypath):
        sub_path = mypath.count(os.path.sep)
        if orig_path + int(options.depth) <= sub_path:
            del dirs[:]
        for subdir in dirs:
            dirname = subdir+"/"
            if (dirname.startswith(options.directory)):
                print(os.path.join(mypath,subdir))
                paths.append(os.path.abspath(os.path.join(mypath,subdir)))
    dot = "."
else: #if (options.directory == ".") or (options.directory == "./"):
    cur_dirn=os.popen('pwd').read()
    cur_dir=cur_dirn.rstrip('\n') 
    paths.append(cur_dir)
    print("here .. no?")
    #time.sleep(1)
    dot = ""


#                                                          CREATING TMP DIRECTORY
rand=int(round(random.random(), 6)*1000000)
tmp="tmpfiles-"+str(rand)
try:
    os.mkdir(tmp)
except:
    print(str(tmp)+" already exists .. weird!")

#                                                          GETTING THE REFERENCE PDB
origin=options.origin
originname=Path(origin).stem
if options.origin != " " :
    if origin.endswith(".pdb"):
        try:
            shutil.copyfile(str(origin), str(originname)+'.pdb')
        except shutil.SameFileError:
            print('reference pdb is already here')
    elif not os.path.isfile(str(originname)+'.pdb'):
        print('============================================\n getting the reference pdb')
        os.system('wget http://www.rcsb.org/pdb/files/'+str(originname)+'.pdb >/dev/null')
#        shutil.move(origin, str(tmp)+'/')


#                                                          GETTING THE PDB
if os.path.isfile(pdb):
    print('============================================\n getting the pdb')
    try:
        shutil.copyfile(str(pdb), str(pdbname)+'.pdb')
    except shutil.SameFileError:
        print('pdb already here')
else:
    print('============================================\n pdb does not exist, will download it')
    try:
        os.system('wget http://www.rcsb.org/pdb/files/'+str(pdbname)+'.pdb >/dev/null')
        os.system('/data/id23eh2/inhouse/opid232/rizk/store/fixpdb.sh '+str(pdbname)+' '+str(LIG))
#        shutil.move(pdb, str(tmp)+'/')
    except:
        print("pdb name is not recognizable .. quitting")
        quit

#                                                          OMIT - APO PDB
if options.not_omit is True:
    notomit = "--full_map"
    maptype = "-full"
    maptype = "full_"
else:
    maptype = ""
    notomit = ""
    maptype = "apo_"
    if os.path.isfile('apo_'+str(pdbname)+'.pdb'):
        print('============================================\n the apo-pdb already exists')
        
#        shutil.copyfile('apo_'+str(pdbname)+'.pdb', str(tmp)+'/apo_'+str(pdbname)+'.pdb')
    else:
        print('============================================\n making the apo of the pdb')
        os.system('/data/id23eh2/inhouse/opid232/rizk/real_data/pyt_getstats/make_apo_lig.py -l '+str(LIG)+' -p '+str(pdb)+'  --apo >/dev/null')

#                                                          MTZ MAP INPUT 
if options.fft_map != "" :
    #shutil.copy(str(options.mtz_map), str(tmp)+'/.')
    mtz_map=os.path.basename(options.fft_map)
    mtz_map_dir=os.path.dirname(options.fft_map)
    mapname=Path(mtz_map).stem
    last_part=mapname.split("_")[-1]
    if last_part == "2fofc":
        map2name=mapname
        mapname=map2name.rstrip(last_part)+"fofc"
        #shutil.copy(mtz_map_dir+'/'+mapname+'.map', str(tmp)+'/.')
    elif last_part == "fofc":
        map2name=mapname.rstrip(last_part)+"2fofc"
        #shutil.copy(mtz_map_dir+'/'+map2name+'.map', str(tmp)+'/.')
    maptype=""
elif options.user_mtz != "":
    mtz_map=os.path.basename(options.user_mtz)
    map2name=Path(mtz_map).stem+'.2fofc'
    mapname=Path(mtz_map).stem+'.fofc'
    maptype=""
else:
    mtz_map=str(options.refine)+"_refined"+str(pdbname)+".mtz"
    pdbRefined=str(options.refine)+"_refined"+str(pdbname)+".pdb"
    map2name=Path(mtz_map).stem+'.2fofc'
    mapname=Path(mtz_map).stem+'.fofc'

#                                                          HKL FILE
if options.hkl is None:
    hkl=""
else:
    hkl="-f "+str(options.hkl)



#                                                        MOVING TO THE TMPFOLDER
os.chdir(tmp)

for folder in sorted(paths):
    print("the folder is : "+str(folder)+"\n")
    #get completeness
    try:
        xscalelp = open(folder+'/XSCALE.LP')
        for line in xscalelp:
            if re.search('total', line):
                c=line.split()[4].strip("%")
                comp="comp_"+c+"_percent"
                break
    except FileNotFoundError:
        print ("Weirdly, XSCALE.LP does not exist ..")
        comp=""
    if (options.fft_map == "") and (options.user_mtz == ""):
        if os.path.isfile(str(folder)+'/'+str(maptype)+str(mtz_map)) and options.overwrite == False : 
            print(' '+str(maptype)+str(mtz_map)+' already exists')
            shutil.copyfile(str(folder)+'/'+str(maptype)+str(mtz_map), str(maptype)+str(mtz_map))
            shutil.copyfile(str(folder)+'/'+str(maptype)+str(pdbRefined), str(maptype)+str(pdbRefined))
        else:
            if os.path.isfile(str(folder)+'/junk_xdsconv.mtz'):
                print('hkl already converted')
                shutil.copyfile(str(folder)+'/junk_xdsconv.mtz', 'junk_xdsconv.mtz')
            else:
                print('============================================\n '+str(maptype)+str(mtz_map)+' doesnt exist\nrunning convert_mtz_hkl.py')
                shutil.copyfile(str(folder)+'/'+str(options.hkl), str(options.hkl))
                os.system('/data/id23eh2/inhouse/opid232/rizk/real_data/pyt_getstats/convert_mtz_hkl.py '+str(hkl)+'>/dev/null')
            print('============================================\n running get_phases_omit-full.py')
            os.system('/data/id23eh2/inhouse/opid232/rizk/real_data/pyt_getstats/get_phases_omit-full_v2.py -l '+str(LIG)+' -p '+str(pdbfullpath)+' --'+str(options.refine)+' '+str(notomit)+' '+str(hkl)+' >/dev/null')
            try:
                shutil.copyfile(str(maptype)+str(mtz_map), str(folder)+'/'+str(maptype)+str(mtz_map))
                shutil.copyfile(str(maptype)+str(pdbRefined), str(folder)+'/'+str(maptype)+str(pdbRefined))
            except FileNotFoundError:
                shutil.copyfile(str(maptype)+str(pdbname)+'_refine_001.mtz', str(folder)+'/'+str(maptype)+str(mtz_map))
                shutil.copyfile(str(maptype)+str(pdbname)+'_refine_001.pdb', str(folder)+'/'+str(maptype)+str(pdbRefined))
    else:
        print (" map or mtz been added by user")
        try:
            shutil.copy(str(folder)+'/'+map2name+'.map', '.')
            shutil.copy(str(folder)+'/'+mapname+'.map', '.')
        except FileNotFoundError:
            shutil.copy(str(folder)+'/'+str(options.user_mtz), '.')
    #if os.path.isfile(str(folder)+'/'+str(maptype)+'refmacout_'+str(pdbname)+'.pdb'):
    #    shutil.copy(str(folder)+'/'+str(maptype)+'refmacout_'+str(pdbname)+'.pdb', '.')

# prepare .pym file

    with open('pym_'+str(LIG)+'.pml', 'w') as outpm:                               
        outpm.write('load ../'+str(pdbname)+'.pdb\n')
    if options.origin != " " :
        print('Origin set by user')
        with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
            outpm.write('load ../'+str(originname)+'.pdb\nsuper '+str(pdbname)+', '+str(originname)+'\ndelete '+str(originname)+'\n')
    else:
        originname=str(maptype)+str(options.refine)+"_refined"+str(pdbname)
        with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
            outpm.write('load '+str(originname)+'.pdb\nsuper '+str(pdbname)+', '+str(originname)+'\ndelete '+str(originname)+'\n')
    with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
        outpm.write('load '+str(maptype)+str(mtz_map)+'\n\
select apomodel, '+str(pdbname)+' and (not resn '+str(LIG)+')\n\
util.cbaw apomodel\n\
select chz, resn '+str(LIG)+'\n\
isomesh map2fofc, '+str(maptype)+map2name+', '+str(options.rmsd2fofc)+', apomodel, carve=1.6\n\
color grey, map2fofc\n\
isomesh map2lig, '+str(maptype)+map2name+', '+str(options.rmsd2fofc)+', chz, carve=1.6\n\
color '+str(options.map_color)+', map2lig\n\
isomesh mapfofc, '+str(maptype)+mapname+', '+str(options.rmsdfofc)+', '+str(pdbname)+'\n\
color tv_green, mapfofc\n\
set mesh_negative_color, red\n\
set mesh_negative_visible, 1\n\
set mesh_negative_visible, 0, map2fofc\n\
set mesh_negative_visible, 0, map2lig\n\
hide cartoon, '+str(pdbname)+'\n\
hide spheres, '+str(pdbname)+'\n\
remove resn hoh\n\
show sticks, '+str(pdbname)+'\n\
set stick_radius, 0.1\n\
show lines, '+str(pdbname)+'\n\
center chz\n\
zoom chz, 2\n\
orient chz\n\
set mesh_width, '+str(options.thick)+'\n\
'+str(viewdata)+' \n')
    if options.color != " " :
        print('\n COLOR set by user')
        with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
            outpm.write('color '+str(options.color)+', chz \n')
    with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
        outpm.write('set ray_shadows=0\nset depth_cue=1\nset ray_trace_fog=1\nset orthoscopic=1\nset antialias=1\nbg_color white\nhide (hydro)\n')
    shutil.copyfile('pym_'+str(LIG)+'.pml', str(folder)+'/pym_'+str(LIG)+'_'+str(maptype)+str(options.refine)+'_'+str(mapname)+'.pml')

# start pymol


    if options.photo is True:
        with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
            outpm.write('png image_'+str(LIG)+', dpi=300')
        print('============================================\n running pymol')
        os.system('pymol -c pym_'+str(LIG)+'.pml >/dev/null ')
        os.system('convert -flatten image_'+str(LIG)+'.png image_'+str(LIG)+'-'+str(os.path.basename(folder))+str(maptype)+str(options.refine)+'.png')
        #        while not os.path.exists('image_'+str(LIG)+'-'+str(os.path.basename(folder))+'.png'):
        #           print('waiting the image..')
        #           time.sleep(5)
        print('============================================\n exporting image')
        if dot == ".":
            #print('============================================\n exporting image')
            shutil.copyfile('image_'+str(LIG)+'-'+str(os.path.basename(folder))+str(maptype)+str(options.refine)+'.png', str(pathlib.Path(folder).parent)+'/image_'+str(LIG)+'-'+str(os.path.basename(folder))+'_'+str(maptype)+str(options.refine)+'_'+str(comp)+'.png')
        else:
            #print('============================================\n exporting image')
            shutil.copyfile('image_'+str(LIG)+'-'+str(os.path.basename(folder))+str(maptype)+str(options.refine)+'.png', str(pathlib.Path(folder))+'/image_'+str(LIG)+'-'+str(os.path.basename(os.path.dirname(os.getcwd())))+'_'+str(maptype)+str(options.refine)+str(comp)+'.png')

    if (options.video is True) or (options.gif is True): 
        multiplier = (10 * 20)/int(options.frames) #duration is 10 seconds here, 20 is the fps
        with open('pym_'+str(LIG)+'.pml', 'a') as outpm:
            outpm.write('mset 1 x'+str(options.frames)+'\nutil.mroll(1,'+str(options.frames)+',1)\nset ray_trace_frames=1\nset cache_frames=0\nmclear\nviewport 640, 600\nmpng snaps_'+str(LIG))
        os.system('pymol -c pym_'+str(LIG)+'.pml >/dev/null')
        for i in range(1, (int(options.frames)+1)):
            n='{num:03d}'.format(num=i)
            os.system('convert -flatten snaps_'+str(LIG)+'0'+str(n)+'.png snaps_'+str(LIG)+'0'+str(n)+'.png')
        os.system('ffmpeg -r 20 -f image2 -s 640x600 -start_number 1 -i snaps_'+str(LIG)+'%04d.png -vframes '+str(options.frames)+' -vcodec libx264 -crf 22  -y -pix_fmt yuv420p video_'+str(LIG)+'.mp4')
        shutil.copyfile('video_'+str(LIG)+'.mp4',  str(folder)+'/video_'+str(LIG)+'.mp4')
        if options.gif is True:
            os.system('ffmpeg -i video_'+str(LIG)+'.mp4 -r 15 -vf "scale=312:-1,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -y gif_'+str(LIG)+'.gif') 
            os.system('ffmpeg -i gif_'+str(LIG)+'.gif -filter:v "setpts='+str(multiplier)+'*PTS" -y '+str(LIG)+'.gif')
            os.remove('video_'+str(LIG)+'.mp4')
            os.remove('gif_'+str(LIG)+'.gif')
            shutil.copyfile(str(LIG)+'.gif',  str(folder)+'/'+str(LIG)+'.gif')
            os.remove(str(LIG)+'.gif')

    #if (os.getcwd().endswith(str(rand))):
     #   print(os.getcwd())
      #  for dele in os.listdir():
       #     os.remove(dele)
            
os.chdir(starting_directory)
#print(os.getcwd())
print("============================================\n removing "+str(tmp))
os.system("rm -r "+tmp)
#print("rm -r "+tmp)
