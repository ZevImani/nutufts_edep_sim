import os,sys
import argparse,signal
import numpy as np
import torch
import gc  # Add garbage collection

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    prog='test_save2petasormdb.py',
    description='Example/Test saving simpledet image into petastorm DB.',
    epilog='Text at the bottom of help')

parser.add_argument('--input-edepsim', '-i', required=True, type=str,
                    help='path to input edepsim file')
parser.add_argument('--pdgcode', '-pdg', required=True, type=int,
                    help='Indicate the particle type being saved')
parser.add_argument("--visualize", '-v', required=False, action='store_true', default=False,
                    help='if flag provided, will visualize the images made before moving to next entry')
parser.add_argument('--start_index', '-s', required=False,type=int,default=0,
                    help='Starting entry index to process from')
parser.add_argument('--batch-size', '-b', required=False,type=int,default=128,
                    help='Number of events per batch (default: 128)')
parser.add_argument('--crop', '-c', required=True, type=int, default=64,
                    help='Integer size of final images to save (default 64x64)')
parser.add_argument('--file_output', '-f', required=False, type=str,
                    help='Folder to save outputs')

# Set arguments and other parameters
args = parser.parse_args()

import numpy as np
from pyspark.sql import SparkSession
from petastorm.etl.dataset_metadata import materialize_dataset
from petastorm.unischema import dict_to_spark_row

# we need ROOT, because EdepSim writes to ROOT by default. So we load it.
import ROOT as rt
rt.gStyle.SetOptStat(0)
# load our library that produces a LArTPC-like image using EDepSim
from simpledet import simpledet
# load schema definition
import simpledet.petastorm.petastorm_schema as schematools

# load the class that executes the conversion
IEdepSim = simpledet.edepsim.EDepSimInterface()
print(IEdepSim)

# input file: output of EDepSim
if not os.path.exists(args.input_edepsim):
    print("Cannot find input EDepSim file at ",args.input_edepsim)
    print("Quitting")
    sys.exit(1)

inputfile = rt.TFile( args.input_edepsim, "open" )
edeptree = inputfile.Get("EDepSimEvents")
if edeptree is None:
    print("Cannot load the expected EDepSimEvents ROOT tree in the input file. Qutting.")
    sys.exit(1)
nentries = edeptree.GetEntries()
print("Loaded EDepSimEvents tree. Number of entries: ",nentries)

cropsize = args.crop

my_images = []
my_images_xz = []
my_moms = []
batch = 0
event_count = 0
bad_events = 0
contains = []

## Starting index and batch size parameters
start_index = args.start_index
batch_size = args.batch_size
print(f"Processing starting from entry {start_index} with batch size {batch_size}")

pdg = args.pdgcode

## number of entries in tree
nentries = edeptree.GetEntries()

# Calculate end index based on start index and how many batches to process
# Process enough entries to fill at least 10 batches from the start index
entries_to_process = batch_size * 10
end_index = min(start_index + entries_to_process, nentries)

if not args.visualize:
    save_folder = args.file_output + "/"
    os.makedirs(save_folder, exist_ok=True)
    if (start_index + entries_to_process) > nentries:
        print("Not enough entries")
        exit( )

print(f"Process entries from {start_index} to {end_index-1} (total: {end_index - start_index} entries)")

testing = False
if testing:
    end_index = nentries

# we collect data for each entry, i.e. simulated image
data = []
for ientry in range(start_index, end_index):

    if ientry > nentries:
        break

    edeptree.GetEntry(ientry)
    print("=========================================")
    print("[ENTRY ",ientry,"]")

    # Get Primary Information
    prim_v = edeptree.Event.Primaries
    print("number of primary vertices: ",prim_v.size())
    prim_mom4_v = []
    prim_pos4_v = []
    prim_pdg_v = []
    for ivertex in range(prim_v.size()):
        print("VERTEX[",ivertex,"]")
        vertex = prim_v.at(ivertex)
        prim_pos = np.zeros(4,dtype=np.float32)
        for v in range(3):
            prim_pos[v] = vertex.GetPosition()[v]*0.1 # mm to cm
        prim_pos[3] = vertex.GetPosition()[3]
        print("  pos4=",prim_pos)
        prim_pos4_v.append( prim_pos )

        mom4_v = []
        pdgcode_v = []
        for iprim in range(vertex.Particles.size()):
            primpart = vertex.Particles.at(iprim)
            prim_mom = np.zeros(4,dtype=np.float32)
            for v in range(4):
                prim_mom[v] = primpart.GetMomentum()[v]
            print("  primary[",iprim,"] ",primpart.GetName()," pdg=",primpart.GetPDGCode()," trackid=",primpart.GetTrackId())
            print("    mom4=",prim_mom)
            pdgcode_v.append( primpart.GetPDGCode() )
            mom4_v.append(prim_mom)
        prim_mom4_v.append( mom4_v )
        prim_pdg_v.append( pdgcode_v )

    if testing:
        mom4 = np.zeros(4,dtype=np.float32)
        mom4 = prim_mom4_v[0][0]
        z,y,x,Etot = mom4
        if 285 < x < 295:
            if 125 < y < 135:
                if 65 < z < 75:
                    print("FOUND YOU")
                    print(z,y,x,Etot)
                else:
                    continue
            else:
                continue
        else:
            continue



    # Get Trajectory information
    traj_v = edeptree.Event.Trajectories
    # print("Number of trajectories: ",traj_v.size())
    #for itraj in range(traj_v.size()):
    #    traj = traj_v.at(itraj)
    #    print("Trajectory[",itraj,"] trackid=",traj.GetTrackId()," nsteps=",traj.Points.size())
    first_trajpoint_pos4 = traj_v.at(0).Points.at(0).GetPosition()
    first_edep = np.zeros(4,dtype=np.float32)
    for v in range(3):
        first_edep[v] = first_trajpoint_pos4[v]*0.1 # mm to cm
    first_edep[3] = first_trajpoint_pos4[3]

    # Get Location of Energy deposits and turn into an image
    seghit_v = edeptree.Event.SegmentDetectors["drift"]
    # print("number of seghits: ",seghit_v.size())

    # some meta data?
    first_seghit = seghit_v.at(0)

    depth = IEdepSim.distance_to_readout_plane # default is 128.0 cm, in the future we can vary this

    # process the edep sim information and make the 2D projection readout plane image
    isgood = IEdepSim.processSegmentHits( seghit_v )

    """
    PyObject* makeNumpyArrayCrop( const TG4HitSegmentContainer& hit_container, int img_pixdim,
                  int offset_x_pixels, int offset_y_pixels, int rand_pix_from_center );
    """

    # make cropped images (bigger, then filter out if event too large)
    threshold = 0.005 # default
    # cropped_dict = IEdepSim.makeNumpyArrayCrop( seghit_v, cropsize, -64, 0, threshold, 0 )
    # cropped_dict = IEdepSim.makeNumpyArrayCrop( seghit_v, cropsize*2, 0, 0, threshold, 0 )

    # cropped_image = cropped_dict["edep"]

    # variables to figure out
    mom4 = np.zeros(4,dtype=np.float32)
    mom4 = prim_mom4_v[0][0]
    # pre_edep_len = np.power( prim_pos4_v[0][:3]-first_edep[:3], 2 ).sum()
    # print("pre_edep_len: ",pre_edep_len)
    # dedx_20pix = np.zeros(20,dtype=np.float32)
    z,y,x,Etot = mom4 ## PILArNet Coordinates

    # entry_data = {"partition":partition_label,
    #               "runid":args.runid,
    #               "entry":ientry,
    #               "pdgcode":prim_pdg_v[0][0],
    #               "depth":depth,
    #               "momentum4":mom4,
    #               "preedeplen":pre_edep_len,
    #               "dedx_20pix":dedx_20pix,
    #               "edepimage":cropped_image}

    # Background threshold
    threshold = 0.005

    # Keep Contained Protons: no pixels outside
    if pdg == 2212:

        # Create extra large image
        big_dict = IEdepSim.makeNumpyArrayCrop( seghit_v, cropsize*2, 0, 0, threshold, 0 )
        big_img = big_dict["edep"]

        H, W = big_img.shape
        center_start = (H - cropsize) // 2
        center_end = center_start + cropsize

        # Check if outer edges have non-zero pixels
        img_edges = np.copy(big_img)
        img_edges[center_start:center_end, center_start:center_end] = 0
        if np.any(img_edges != 0):
            bad_events += 1

            # Clean up memory for this iteration before continuing
            del big_dict, big_img, img_edges
            gc.collect()
            continue

        # Keep center crop of image
        img = big_img[center_start:center_end, center_start:center_end]

        # xz projection via C++ pipeline (same diffusion model as xy)
        IEdepSim.processSegmentHitsXZ( seghit_v )
        big_dict_xz = IEdepSim.makeNumpyArrayCrop( seghit_v, cropsize*2, 0, 0, threshold, 0 )
        big_img_xz = big_dict_xz["edep"]
        img_xz = big_img_xz[center_start:center_end, center_start:center_end]

        if testing:
            plt.imshow(img, cmap='gray')
            plt.title(mom4)
            plt.savefig("zzz_2.png")
            exit()

        # Clean up intermediate variables
        del big_dict, big_img, img_edges, big_dict_xz, big_img_xz

    # Keep Contained Electrons: 90% energy contained in image
    if pdg == 11:

        # full detector
        full_dict = IEdepSim.makeNumpyArrayCrop( seghit_v, 512, 0, 0, threshold, 0 )
        full_img = full_dict["edep"]

        # final image
        # img_dict = IEdepSim.makeNumpyArrayCrop( seghit_v, cropsize, 0, 0, threshold, 0 )
        # img = img_dict["edep"]

        H, W = full_img.shape
        center_start = (H - cropsize) // 2
        center_end = center_start + cropsize
        img = full_img[center_start:center_end, center_start:center_end]

        # Determine energies
        e_tot = np.sum(full_img)
        e_img = np.sum(img)

        # percentage of contained energy
        e_contained = e_img / e_tot
        if e_tot == 0:
            e_contained = 0

        # skip events with over exposure (TODO: fix this)
        if np.isnan(e_contained):
            bad_events += 1

            # Clean up memory for this iteration before continuing
            del full_dict, full_img
            gc.collect()
            continue

        # only keep 80% contained events
        if e_contained < 0.8:
            bad_events += 1

            # Clean up memory for this iteration before continuing
            del full_dict, full_img
            gc.collect()
            # continue

        contains.append(e_contained)

        # xz projection via C++ pipeline (same diffusion model as xy)
        IEdepSim.processSegmentHitsXZ( seghit_v )
        full_dict_xz = IEdepSim.makeNumpyArrayCrop( seghit_v, 512, 0, 0, threshold, 0 )
        full_img_xz = full_dict_xz["edep"]
        img_xz = full_img_xz[center_start:center_end, center_start:center_end]

        # Clean up intermediate variables
        try:
            del full_dict, full_img, full_dict_xz, full_img_xz
        except:
            pass

    ## Log image and momentum
    my_images.append(img)
    my_images_xz.append(img_xz)
    my_moms.append(np.array([x,y,z]))
    event_count += 1

    ## Plot grid of events
    if args.visualize and len(my_images) >= 16:

        grid_size = 4
        fig, axes = plt.subplots(grid_size, grid_size, figsize=(8, 8))
        axes = axes.ravel() # Flatten axes array for easy iteration

        print(len(my_images))

        # fig.suptitle("Edep-Sim: p="+str(np.round(my_moms[0],1)), fontsize=20)
        fig.suptitle("Edep-Sim", fontsize=20)
        for i in range(len(axes)):
            # try:

            this_img = my_images[i]
            # this_img[this_img != 0] = 10 # highlight nonzeros
            axes[i].imshow(this_img , cmap='gray')
            axes[i].axis('off')
            # axes[i].set_title(str(np.round(my_moms[i],1)))
            if len(contains) > 0:
                axes[i].set_title(str(np.round(contains[i],3)))
            else:
                axes[i].set_title(str(np.round(my_moms[i],2)))

            # except:
            #     break
            # if i in [1,4,7]:
            #     axes[i].set_title(str(np.round(my_moms[i],2)))

        plt.tight_layout()
        plt.savefig("zzz.png")
        print("Saved: zzz.png")
        print("Image size:", my_images[i].shape)
        print("Bad events", bad_events)

        if len(contains) > 0:
            contains = np.array(contains)
            print(np.min(contains), np.mean(contains), np.max(contains))

        # Clean up matplotlib objects
        plt.close(fig)
        del fig, axes

        exit()

    ## Save batch
    else:

        if event_count == batch_size:

            print(f"Saving batch {batch} with {len(my_images)} events")

            # Calculate batch filename using start_index for uniqueness
            batch_filename = f"batch_{start_index//batch_size + batch}.npy"
            xz_filename   = f"batch_xz_{start_index//batch_size + batch}.npy"
            mom_filename   = f"batch_mom_{start_index//batch_size + batch}.npy"

            np.save(f"{save_folder}/{batch_filename}", my_images)
            np.save(f"{save_folder}/{xz_filename}", my_images_xz)
            np.save(f"{save_folder}/{mom_filename}", my_moms)
            print(f"Saved: {batch_filename}, {xz_filename}, and {mom_filename}")

            batch += 1

            # Clear batch data to free memory
            my_moms = []
            my_images = []
            my_images_xz = []
            event_count = 0

    # MEMORY CLEANUP AFTER EACH ITERATION
    # Clear local variables that might hold references to ROOT objects
    try:
        del prim_v, prim_mom4_v, prim_pos4_v, prim_pdg_v
        del traj_v, first_trajpoint_pos4, first_edep
        del seghit_v, first_seghit
        del mom4, dedx_20pix
        if 'img' in locals():
            del img
        if 'img_xz' in locals():
            del img_xz
        if 'vertex' in locals():
            del vertex
        if 'mom4_v' in locals():
            del mom4_v
        if 'pdgcode_v' in locals():
            del pdgcode_v
        if 'primpart' in locals():
            del primpart
        if 'prim_mom' in locals():
            del prim_mom
        if 'prim_pos' in locals():
            del prim_pos
    except:
        pass

    # Force garbage collection every iteration to prevent memory accumulation
    gc.collect()

    # Additional ROOT memory cleanup every 10 iterations
    if ientry % 10 == 0:
        rt.gROOT.GetListOfCanvases().Delete()
        rt.gROOT.GetListOfBrowsers().Delete()
        rt.gSystem.ProcessEvents()

    # Print memory usage periodically
    if ientry % 50 == 0:
        print(f"Processed {ientry} entries, current batch: {batch}, bad events: {bad_events}")

# Save any remaining events in the final partial batch
# if not args.visualize and event_count > 0:

#     save_folder = args.file_output + "/"

#     os.makedirs(save_folder, exist_ok=True)

#     print(f"Saving to {save_folder}")
#     print(f"Saving final partial batch {batch} with {len(my_images)} events")
#     batch_filename = f"batch_{start_index//batch_size + batch}.npy"
#     mom_filename = f"batch_mom_{start_index//batch_size + batch}.npy"

    # np.save(f"{save_folder}/{batch_filename}", my_images)
    # np.save(f"{save_folder}/{mom_filename}", my_moms)
#     print(f"Saved final batch: {batch_filename} and {mom_filename}")

# Final cleanup before exiting
print("DONE")
print("Bad events", bad_events)
print("Containment:", contains)
if len(contains) > 0:
    contains = np.array(contains)
    print(np.min(contains), np.mean(contains), np.max(contains))

# Final memory cleanup
inputfile.Close()
del inputfile, edeptree
gc.collect()

# Clear ROOT global objects
rt.gROOT.GetListOfCanvases().Delete()
rt.gROOT.GetListOfBrowsers().Delete()
rt.gROOT.GetListOfFiles().Delete()
