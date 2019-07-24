


use arcstar::sae_types::*;
//use rayon::prelude::*;

//use image::{ RgbImage };
//use imageproc::drawing;

use graphics_buffer::*;
//use graphics::types::{Color, ColorComponent};

//use rand::{Rng};

use crate::slab::{SlabStore};
//use num_traits::float::Float;

pub struct FeatureTracker {
    store: SlabStore,
}

impl FeatureTracker {
    pub const RAINBOW_WHEEL_DIM: usize = 12;
    pub const RAINBOW_WHEEL_GEN:[ [u8; 3]; Self::RAINBOW_WHEEL_DIM] = [
            [255, 0, 0 ],
            [255, 127, 0],
            [255, 255, 0],
            [127, 255, 0],
            [0, 255, 0],
            [0, 255, 127],
            [0, 255, 255],
            [0, 127, 255],
            [0, 0, 255],
            [127, 0, 255],
            [255, 0, 255],
            [255, 0, 127],
            ];

    pub fn new() -> FeatureTracker {
        FeatureTracker {
            store:  SlabStore::new(),
        }
    }

    pub fn add_and_match_feature(&mut self, new_evt: &SaeEvent) -> Option<SaeEvent> {
        let res = self.store.add_and_match_feature(new_evt);
        match res {
            Some(feat) => Some(   feat.event),
            None => None
        }
    }


    fn render_track(segments:&Vec<[f64;4]>, buffer: &mut RenderBuffer, color_hint: usize) {
      for seg in segments {
        Self::render_one_line(*seg, buffer, color_hint);
      }
    }

    fn render_one_line(line: [f64; 4], buffer: &mut RenderBuffer, color_hint: usize)  {
      let color_idx =  (color_hint * (FeatureTracker::RAINBOW_WHEEL_DIM / 3) + color_hint ) % FeatureTracker::RAINBOW_WHEEL_DIM;
      let rgb_data = FeatureTracker::RAINBOW_WHEEL_GEN[color_idx];
      let px:[f32; 4] = [
        (rgb_data[0] as f32)/255.0,
        (rgb_data[1] as f32)/255.0,
        (rgb_data[2] as f32)/255.0,
        1.0 //alpha
      ];
      graphics::line(
        px,
        1.0,
        line,
        IDENTITY,
        buffer
      );

    }



    pub fn render_tracks_to_file(&self, nrows: u32, ncols: u32, lead_time_horizon: SaeTime, out_path: &str )  {
      let mut buffer = RenderBuffer::new(ncols, nrows);
      buffer.clear([0.0, 0.0, 0.0, 1.0]);

      let chains_list:Vec<Vec<[f64; 4]>> = (0..nrows*ncols).fold(vec!(), |mut acc, idx| {
        let row: u16 = (idx / ncols) as u16;
        let col: u16 = (idx % ncols) as u16;
        let chain = self.store.chain_for_point(row, col, 0);
        if chain.len() > 2 { //TODO arbitrary cutoff
          let lead_evt = &chain[0];
          let mut total_dist_x = 0;
          let mut total_dist_y = 0;

          if lead_evt.timestamp >= lead_time_horizon {
            //add all events in the chain
            let mut chain_vec:Vec<[f64; 4]> = vec!();
            for i in 0..(chain.len()-1) {
              let evt = &chain[i];
              let old_evt = &chain[i+1];
              let dist_x = (evt.col as i32) - (old_evt.col as i32);
              let dist_y = (evt.row as i32) - (old_evt.row as i32);
              total_dist_x += dist_x;
              total_dist_y += dist_y;


              let le_line: [f64; 4] = [evt.col as f64, evt.row as f64, old_evt.col as f64, old_evt.row as f64];
              chain_vec.push(le_line);
            }
            acc.push(chain_vec);

            //if 0 == total_dist_y {
            //  println!("total_dist xy: {}, {}", total_dist_x, total_dist_y);
            //}
          }
        }
        acc
      });

      println!("rendering {} lines", chains_list.len());
      for (i, chain) in chains_list.iter().enumerate() {
        Self::render_track(chain, &mut buffer, i);
      }


      buffer.save(out_path).expect("Couldn't save");

    }

}

//
//#[cfg(test)]
//mod tests {
//    use super::*;
//x
//    #[test]
//    fn test_find_neighbors() {
//        let mut tracker = FeatureTracker::new();
//        //test empty tracker
//    }
//}