#!/bin/bash

# convert multiple png files to a video file
# requires that libx264 and ffmpeg are installed
rm sae_surface.mp4

ffmpeg -r 24 \
 -f image2 \
 -start_number 1 \
 -i ./out/sae_surf_%04d.png \
 -vcodec libx264 \
 -crf 25 \
 -pix_fmt yuv420p \
 sae_surface.mp4
