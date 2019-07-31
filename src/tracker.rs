
use arcstar::sae_types::*;
use arcstar::detector;
use crate::slab;

use image::{ RgbImage };
use imageproc::drawing;

use rand::Rng;


pub struct FeatureTracker {
  n_pixel_rows: u32,
  n_pixel_cols: u32,
  store: slab::SlabStore,
  sae_rise: SaeMatrix,
  sae_fall: SaeMatrix,
  ///related to tmax from ACE paper
  time_window: SaeTime,
}

impl FeatureTracker {
  pub const RAINBOW_WHEEL_DIM: usize = 12;
  pub const RAINBOW_WHEEL_GEN: [[u8; 3]; Self::RAINBOW_WHEEL_DIM] = [
    [255, 0, 0],
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

  pub fn new(img_w: u32, img_h: u32, time_window: SaeTime) -> FeatureTracker {
    FeatureTracker {
      n_pixel_cols: img_w,
      n_pixel_rows: img_h,
      store: slab::SlabStore::new(),
      sae_rise: SaeMatrix::zeros(img_h as usize, img_w as usize),
      sae_fall: SaeMatrix::zeros(img_h as usize, img_w as usize),
      time_window: time_window,
    }
  }


  /// process a single pixel change event,
  /// return it as a feature if it's a corner
  pub fn process_one_event(&mut self, evt: &SaeEvent) -> Option<SaeEvent> {
    let row = evt.row as usize;
    let col = evt.col  as usize;

    let sae_pol =  match evt.polarity {
      1 => &mut self.sae_rise,
      0 | _ => &mut self.sae_fall,
    };

    sae_pol[(row, col)] = evt.timestamp;
    let feature =  detector::detect_and_compute_one(sae_pol, evt);

    if feature.is_some() {
      // this event is a feature -- add it to the store
      let new_evt = feature.clone().unwrap();
      self.track_feature(&new_evt);
    }

    feature
  }

  /// process a list of pixel change events.  these may take place over any time span
  pub fn process_events(&mut self, event_list: &Vec<SaeEvent>) -> Vec<SaeEvent>  {
    let corners: Vec<SaeEvent> = event_list.iter().filter_map(|evt| {
      self.process_one_event(evt)
    }).collect();

    corners
  }

  /// add the given feature to the store, and attempt to match it against other stored features
  fn track_feature(&mut self, new_evt: &SaeEvent) -> Option<SaeEvent> {
    let time_horizon: SaeTime = new_evt.timestamp.max(self.time_window) - self.time_window;
    self.store.add_and_match_feature(new_evt, time_horizon)
  }

  /// Render a representation of the SAE to a file
  pub fn render_sae_frame_to_file(&self,  time_horizon: SaeTime, out_path: &str ) {
    let out_img = self.render_sae_frame(time_horizon );
    out_img.save(out_path).expect("Couldn't save");
  }

  /// Render a representation of the SAE with the given cutoff time horizon
  pub fn render_sae_frame(&self,  time_horizon: SaeTime) -> RgbImage {
    let mut out_img = RgbImage::new(self.n_pixel_cols, self.n_pixel_rows);

    for row in 0..self.n_pixel_rows {
      for col in 0..self.n_pixel_cols {
        let sae_rise_val: SaeTime = self.sae_rise[(row as usize, col as usize)];
        let sae_fall_val: SaeTime = self.sae_fall[(row as usize, col as usize)];
        let sae_val = sae_rise_val.max(sae_fall_val);

        if 0 != sae_val && sae_val >= time_horizon {
          //let total_val = (sae_fall_val + sae_rise_val) as f32;
          //let red_val: u8 = (255.0 * (sae_rise_val as f32) / total_val) as u8;
          //let green_val: u8 = (255.0  * (sae_fall_val as f32) / total_val) as u8;
          let lum:u8 = 255; //(255.0 * sae_pix_frac ) as u8;
          let red_val:u8 =  if sae_rise_val > sae_fall_val { lum } else {0};
          let green_val: u8 = if sae_fall_val > sae_rise_val  { lum } else {0};
          let px_data: [u8; 3] = [red_val, green_val, 0];
          let px = image::Rgb(px_data);

          out_img.put_pixel(col as u32, row as u32, px);
        }

      }
    }


    out_img
  }

  pub const RED_PIXEL: [u8; 3] = [255u8, 0, 0];
  pub const GREEN_PIXEL: [u8; 3] = [0, 255u8, 0];
  pub const YELLOW_PIXEL: [u8; 3] = [255u8, 255u8, 0];
  pub const BLUE_PIXEL: [u8; 3] = [0,0,  255u8];


  /// render all tracks into an image
  pub fn render_tracks(&self, time_horizon: SaeTime) -> RgbImage {
    let nrows = self.n_pixel_rows;
    let ncols = self.n_pixel_cols;
    let mut out_img =   RgbImage::new(ncols , nrows);

    // collect the list of tracks to render
    let chains_list: Vec<Vec<((f32, f32), (f32, f32))>> =
      (0..nrows * ncols).fold(vec!(), |mut acc, idx| {
        let row: u16 = (idx / ncols) as u16;
        let col: u16 = (idx % ncols) as u16;
        let chain = self.store.chain_for_point(row, col, time_horizon);
        if chain.len() > 4 { //TODO arbitrary cutoff
          //add all events in the chain that are after the time horizon
          let mut chain_vec: Vec<((f32, f32), (f32, f32))> = Vec::with_capacity(chain.len());
          for i in 0..(chain.len() - 1) {
            let evt = &chain[i];
            let old_evt = &chain[i + 1];
            let pair = ((evt.col as f32, evt.row as f32), (old_evt.col as f32, old_evt.row as f32));
            chain_vec.push(pair);
          }
          acc.push(chain_vec);
        }
        acc
      });


    //draw each chain in its own color
    for (i, chain) in chains_list.iter().enumerate() {
      let rgb_data = FeatureTracker::RAINBOW_WHEEL_GEN[(i as usize) % FeatureTracker::RAINBOW_WHEEL_DIM];
      let px = image::Rgb(rgb_data);
      for segment in chain {
        drawing::draw_line_segment_mut(&mut out_img, segment.0, segment.1, px);
      }
    }

    out_img
  }

  /// Render lines for all the valid tracks to the given file path
  pub fn render_tracks_to_file(&self, time_horizon: SaeTime, out_path: &str) {
    let out_img = self.render_tracks(time_horizon );
    out_img.save(out_path).expect("Couldn't save");
  }


  /// render events into an image result
  pub fn render_events(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3] ) -> RgbImage {
    let mut out_img =   RgbImage::new(self.n_pixel_cols , self.n_pixel_rows);

    for evt in events {
      let px = match evt.polarity {
        1 => image::Rgb(*rising_pix),
        0 => image::Rgb(*falling_pix),
        _ => unreachable!()
      };

      out_img.put_pixel(evt.col as u32, evt.row as u32, px);
    }

    out_img
  }

  /// render events into an image file
  pub fn render_events_to_file(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3], out_path: &str ) {
    let out_img = self.render_events(events, rising_pix, falling_pix);
    out_img.save(out_path).expect("Couldn't save");
  }

  /// render events as corners to an image result
  pub fn render_corners(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3] ) -> RgbImage {
    let mut out_img =   RgbImage::new(self.n_pixel_cols , self.n_pixel_rows);

    for evt in events {
      let px = match evt.polarity {
        1 => image::Rgb(*rising_pix),
        0 => image::Rgb(*falling_pix),
        _ => unreachable!()
      };

      drawing::draw_cross_mut(&mut out_img, px, evt.col as i32, evt.row as i32);
    }

    out_img
  }

  /// render events as corners to a file
  pub fn render_corners_to_file(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3], out_path: &str ) {
    let out_img = self.render_corners(events, rising_pix, falling_pix);
    out_img.save(out_path).expect("Couldn't save");

  }

//  squiggle that triggers the corner detector, but not clear what this corresponds to IRL
//  const SAMPLE_CORNER_DIM: usize = 8;
//  const SAMPLE_CORNER_GEN: [[i32; 2]; Self::SAMPLE_CORNER_DIM] = [
//    [1, -4], [2, -3], [3, -2], [4, -1],
//    [1, -3], [2, -2], [3, -1],
//    [0, 0]
//  ];

//  //true corner L
//  const SAMPLE_CORNER_DIM: usize = 10;
//  const SAMPLE_CORNER_GEN: [[i32; 2]; Self::SAMPLE_CORNER_DIM] = [
//    [4, 0], [3,0], [2,0], [1,0],
//    [1,-1],
//    [0, -4], [0, -3], [0, -2], [0, -1],
//    [0, 0]
//  ];

  //glider
const SAMPLE_CORNER_DIM: usize = 7;
const SAMPLE_CORNER_GEN: [[i32; 2]; Self::SAMPLE_CORNER_DIM] = [
  [3, 3], [3, -3],
  [2, 2], [2 , -2],
  [1, 1], [1, -1],
  [0, 0],

];



  const SYNTH_TIME_INCR:SaeTime = 10;
  /// Insert a synthetic event at the given position:
  /// append the synthesized event to the modified event list
  pub fn insert_one_synth_event(ctr_row: i32, ctr_col: i32, timestamp: &mut SaeTime, polarity: u8, event_list: &mut Vec<SaeEvent> ) {
    for j in 0..Self::SAMPLE_CORNER_DIM {
      *timestamp += Self::SYNTH_TIME_INCR;
      let dxy = Self::SAMPLE_CORNER_GEN[j];
      let evt_row = ctr_row + dxy[1];
      let evt_col = ctr_col + dxy[0];

      let mut rng = rand::thread_rng();
      let mut ndesc:NormDescriptor = [0.0; NORM_DESCRIPTOR_LEN];

      for i in 0..ndesc.len() {
        ndesc[i] = rng.gen::<f32>();
      }

      let evt = SaeEvent {
        row: evt_row as u16,
        col: evt_col as u16,
        polarity: polarity,
        timestamp: *timestamp,
        norm_descriptor: Some(Box::new(ndesc)),
      };

      event_list.push(evt);
    }
  }


  /// Generate a series of synthetic events for stimulating the tracker
  pub fn process_synthetic_events(img_w: u32, img_h: u32, render_out: bool) {
    const TMAX_FORGETTING_TIME: SaeTime = (0.1 / 1E-6) as SaeTime;
    let mut tracker = Box::new(FeatureTracker::new(img_w, img_h, TMAX_FORGETTING_TIME));
    let mut timestamp: SaeTime = 1;

    let spacer = (2 * slab::MAX_RL_SPATIAL_DISTANCE) as i32;
    let mut ctr_col: i32 = (img_w as i32) - spacer;
    let ctr_row: i32 = (img_w / 2) as i32;
    let mut chunk_count = 0;

    //these three constants in essence define the minimum distinguishable corner grid
    const COL_DECR: i32 = 2; //how many pixels decrement per iteration
    let row_gap: i32 = spacer;
    let col_gap: i32 = spacer;

    let half_width: i32 = (img_w / 2) as i32;
    let total_frames = half_width / COL_DECR;

    let max_time_delta: SaeTime = 24 * (Self::SYNTH_TIME_INCR * Self::SAMPLE_CORNER_DIM as u32);

    for _i in 0..total_frames {
      let mut event_list: Vec<SaeEvent> = vec!();
      chunk_count += 1;

      Self::insert_one_synth_event(ctr_row - row_gap, ctr_col, &mut timestamp, 1,&mut event_list);
      Self::insert_one_synth_event(ctr_row, ctr_col, &mut timestamp, 1, &mut event_list);
      Self::insert_one_synth_event(ctr_row + row_gap, ctr_col, &mut timestamp, 1, &mut event_list);

      Self::insert_one_synth_event(ctr_row - row_gap, ctr_col - col_gap, &mut timestamp, 0, &mut event_list);
      Self::insert_one_synth_event(ctr_row, ctr_col - col_gap, &mut timestamp, 0, &mut event_list);
      Self::insert_one_synth_event(ctr_row + row_gap, ctr_col - col_gap, &mut timestamp,  0,&mut event_list);


      let corners = tracker.process_events( &event_list);
      println!("chunk {} events {} corners {}", chunk_count, event_list.len(), corners.len());

      if render_out {
        let render_corners = true;
        let render_tracks = true;
        let horizon = timestamp.max( max_time_delta) - max_time_delta;

        if render_corners {
          let out_path = format!("./out/sae_{:04}_corners.png", chunk_count);
          tracker.render_corners_to_file( &corners, &FeatureTracker::YELLOW_PIXEL, &FeatureTracker::GREEN_PIXEL, &out_path );
        }
        if render_tracks {
          let out_path= format!("./out/sae_{:04}_tracks.png", chunk_count);
          tracker.render_tracks_to_file(horizon, &out_path);
        }

        let out_path=  format!("./out/saesurf_{:04}.png", chunk_count);
        tracker.render_sae_frame_to_file(horizon, &out_path);

      }

      //assert_eq!(corners.len(), 6 );

      ctr_col -= COL_DECR;
    }
  }


}



//#[cfg(test)]
//mod tests {
//  use super::*;
//  use std::fs::create_dir_all;
//  use std::path::Path;
//
//
//  #[test]
//  fn test_synthetic_worm_race() {
//    let img_w = 320;
//    let img_h = 320;
//    create_dir_all(Path::new("./out/")).expect("Couldn't create output dir");
//    FeatureTracker::process_synthetic_events(img_w, img_h, true);
//  }
//}