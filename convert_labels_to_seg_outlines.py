import skimage.segmentation
import sys
import tifffile

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} LABEL_MASK_IN.ome.tif SEG_OUTLINES_OUT.ome.tif")
    sys.exit(1)

labels_in = sys.argv[1]
seg_out = sys.argv[2]

labels = tifffile.imread(labels_in)
seg = skimage.segmentation.find_boundaries(labels)
seg = skimage.img_as_ubyte(seg)
tifffile.imwrite(seg_out, seg, compression='deflate', tile=(1024, 1024))
