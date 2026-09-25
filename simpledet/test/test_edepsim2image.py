import os,sys
import argparse

parser = argparse.ArgumentParser(
    prog='test_edepsim2image',
    description='Example/Test for turning EdepSim output into a simpledet image',
    epilog='Text at the bottom of help')

parser.add_argument('--input-edepsim', '-i', required=True, type=str,
                    help='path to input edepsim file')
parser.add_argument('--output-format', '-f', required=True, type=str,
                    help='output format. options: [root,petastorm]')
parser.add_argument('--output-rootfile',type=str,
                    help='If output-format is "root", this is required to set the output file')
parser.add_argument('--petastorm-db-folder', '-db', type=str,
                    help='If output-format is "petastorm", this argument is required in order to set the location to the petastorm database folder.')
parser.add_argument('--root-batchmode', '-b', type=bool, default=False, #action=store_true,
                    help='If using "root" mode, run without drawing to canvas.')


args = parser.parse_args()

# we need ROOT, because EdepSim writes to ROOT by default. So we load it.
import ROOT as rt
# load our library that produces a LArTPC-like image using EDepSim
from simpledet import simpledet

# load the class that executes the conversion
IEdepSim = simpledet.edepsim.EDepSimInterface()
print(IEdepSim)

# input file: output of EDepSim
if not os.path.exists(args.input_edepsim):
    print("Cannot find input EDepSim file at ",args.input_edepsim)
    print("Quitting")
    exit()

inputfile = rt.TFile( args.input_edepsim, "open" )
edeptree = inputfile.Get("EDepSimEvents")
if edeptree is None:
    print("Cannot load the expected EDepSimEvents ROOT tree in the input file. Qutting.")
    exit()


nentries = edeptree.GetEntries()
print("Loaded EDepSimEvents tree. Number of entries: ",nentries)

# Setup the output mode.
# We will try ROOT's PyTorch Batch Loader ...
# https://root.cern/doc/v630/RBatchGenerator__PyTorch_8py.html
# we can either dump to ROOT file, or store into a petastorm (local) database
# if args.output_format=="root":
#     if True:
#         print("Preferred way for now is petastorm loader. It's known to work with multithreading, so is expected to be fast.")
#         sys.exit(1)
#     if os.path.exists( args.output_rootfile ):
#         print("output file already exists. will not overwite. output path given:",args.output_rootfile )
#         exit()

#     out_rootfile = rt.TFile( args.output_rootfile, 'new' )
# elif args.output_format=="petastorm":
#     pass
# else:
#     print("Unrecgonized output format")
#     sys.exit(1)


#c1 = rt.TCanvas("c1","",800,600)

for ientry in range(nentries):

    edeptree.GetEntry(ientry)

    seghit_v = edeptree.Event.SegmentDetectors["drift"]
    print("number of seghits: ",seghit_v.size())
    himage = IEdepSim.processSegmentHits( seghit_v )
    print("himage: ",himage)
    himage.Draw("colz")
    c1.Update()
    print("[ENTER] to continue")
    input()





