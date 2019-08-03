
use arcstar::sae_types::*;
use arcstar::detector;
use crate::slab;

use image::{ RgbImage };
use imageproc::drawing;



struct ChainContainer {
  /// the line segments to render
  chain: Vec<((f32, f32), (f32, f32))>,
  /// color for this chain
  color: [u8; 3],
}


///
/// Note that portions of this work are based on
/// "Asynchronous Corner Detection and Tracking for Event Cameras in Real-Time"
/// by Alzugaray et al, which we'll refer to as the "ACD paper" herein,
///
/// and "ACE: An Efficient Asynchronous Corner Tracker for Event Cameras" by Alzugaray et al,
/// which we'll refer to as the "ACE paper"
///
///
pub struct FeatureTracker {
  /// The number of rows of pixels provided by our sensor
  n_pixel_rows: u32,
  /// The number of columns of pixels provided by our sensor
  n_pixel_cols: u32,
  /// Data store for tracking detected features
  store: slab::SlabStore,
  /// Surface of Active Events (SAE) for pixel changes rising about a threshold
  sae_rise: SaeMatrix,
  /// Surface of Active Events (SAE) for pixel changes falling below a threshold
  sae_fall: SaeMatrix,
  /// Related to ∆tmax from ACE paper -- the sliding time window of event relevance
  time_window: SaeTime,
  /// Related to (κ = 50ms) event redundancy filter from ACD paper
  ref_time_filter: SaeTime,
}

impl FeatureTracker {

  /// Some convenient pixel color constants for rendering tracker progress
  pub const RED_PIXEL: [u8; 3] = [255u8, 0, 0];
  pub const GREEN_PIXEL: [u8; 3] = [0, 255u8, 0];
  pub const YELLOW_PIXEL: [u8; 3] = [255u8, 255u8, 0];
  pub const BLUE_PIXEL: [u8; 3] = [0, 0, 255u8];
  pub const MAGENTA_PIXEL: [u8; 3] = [255u8, 0, 255u8];

  pub fn new(img_w: u32, img_h: u32, time_window: SaeTime, ref_time_filter: SaeTime) -> FeatureTracker {
    FeatureTracker {
      n_pixel_cols: img_w,
      n_pixel_rows: img_h,
      store: slab::SlabStore::new(),
      sae_rise: SaeMatrix::zeros(img_h as usize, img_w as usize),
      sae_fall: SaeMatrix::zeros(img_h as usize, img_w as usize),
      time_window: time_window,
      ref_time_filter: ref_time_filter,
    }
  }

  /// Returns a possibly empty list of events that were detected as corners.
  ///
  /// # Arguments
  ///
  /// * `event_list` - A list of pixel change events. These may occur over any time span.
  ///
  pub fn process_events(&mut self, event_list: &Vec<SaeEvent>) -> Vec<SaeEvent>  {
    let corners: Vec<SaeEvent> = event_list.iter().filter_map(|evt| {
      self.process_one_event(evt)
    }).collect();

    corners
  }

  /// Returns a modified event with descriptor if the event was detected as a corner.
  ///
  /// # Arguments
  ///
  /// * `event` - A single pixel change event.
  ///
  pub fn process_one_event(&mut self, event: &SaeEvent) -> Option<SaeEvent> {
    let row = event.row as usize;
    let col = event.col  as usize;

    let sae_pol =  match event.polarity {
      1 => &mut self.sae_rise,
      0 | _ => &mut self.sae_fall,
    };

    //filter noisy/redundant events using reference time threshold
    let last_timestamp = sae_pol[(row, col)];
    if (event.timestamp - last_timestamp) < self.ref_time_filter {
      //ignore redundant events (noisy event filter)
      return None
    }

    sae_pol[(row, col)] = event.timestamp;
    let feature =  detector::detect_and_compute_one(sae_pol, event);

    if feature.is_some() {
      // this event is a feature -- add it to the store
      let new_feature = feature.clone().unwrap();
      self.track_feature(&new_feature);
    }

    feature
  }


  /// Returns a matched parent feature for the given feature, if a match exists
  ///
  /// Add the given feature to the store, and attempt to match it against other stored features
  fn track_feature(&mut self, new_evt: &SaeEvent) -> Option<SaeEvent> {
    let time_horizon: SaeTime = new_evt.timestamp.max(self.time_window) - self.time_window;
    self.store.add_and_match_feature(new_evt, time_horizon)
  }


  /// Returns a list of valid tracks
  ///
  /// # Arguments
  ///
  /// * `track_len_filter` - Tracks must be at least this long (minimum of 2 enforced)
  ///
  fn collect_tracks(&self, track_len_filter: usize, time_horizon: SaeTime) -> Vec<ChainContainer> {
    let nrows = self.n_pixel_rows;
    let ncols = self.n_pixel_cols;
    //it takes at least two features to make a track:
    // enforce track_len_filter > 1
    let track_len_filter = if track_len_filter > 1 { track_len_filter } else { 2 };

    // collect the list of tracks to render
    let chains_list: Vec<ChainContainer> =
      (0..nrows * ncols).fold(vec!(), |mut acc, idx| {
        let row: u16 = (idx / ncols) as u16;
        let col: u16 = (idx % ncols) as u16;
        if let Some(track) = self.store.track_for_point(row, col, time_horizon) {
          let track_len = track.chain.len();
          if track_len >= track_len_filter {
            //add all events in the chain that are after the time horizon
            let mut chain_vec: Vec<((f32, f32), (f32, f32))> = Vec::with_capacity(track_len);

            for i in 0..(track_len - 1) {
              let evt = &track.chain[i];
              let old_evt = &track.chain[i + 1];
              let pair = ((evt.col as f32, evt.row as f32), (old_evt.col as f32, old_evt.row as f32));
              chain_vec.push(pair);
            }

            acc.push(ChainContainer {
              chain: chain_vec,
              color: track.color,
            });
          }
        }
        acc
      });

    chains_list
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
          let lum:u8 = 255; //(255.0 * sae_pix_frac ) as u8;
          let red_val:u8 =  if sae_rise_val > sae_fall_val { lum } else { 0 };
          let green_val: u8 = if sae_fall_val > sae_rise_val  { lum } else { 0 };
          let px_data: [u8; 3] = [red_val, green_val, 0];
          let px = image::Rgb(px_data);

          out_img.put_pixel(col as u32, row as u32, px);
        }

      }
    }

    out_img
  }

  /// Render all tracks into an image
  pub fn render_tracks(&self, time_horizon: SaeTime) -> RgbImage {
    let chains_list = self.collect_tracks(3, time_horizon);
    self.render_track_segments(chains_list)
  }

  fn render_track_segments(&self, chains_list: Vec<ChainContainer>) -> RgbImage {
    let nrows = self.n_pixel_rows;
    let ncols = self.n_pixel_cols;
    let mut out_img =   RgbImage::new(ncols , nrows);

    //draw each chain in its own color
    for chain_container in chains_list.iter() {
      let px = image::Rgb(chain_container.color);

      for segment in chain_container.chain.iter() {
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

  /// Render events into an image result
  pub fn render_events(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3] ) -> RgbImage {
    let mut out_img =   RgbImage::new(self.n_pixel_cols , self.n_pixel_rows);

    for evt in events {
      let px = match evt.polarity {
        1 => image::Rgb(*rising_pix),
        0 | _ => image::Rgb(*falling_pix),
      };

      // put a single pixel at every event location
      out_img.put_pixel(evt.col as u32, evt.row as u32, px);
    }

    out_img
  }

  /// Render events into an image file
  pub fn render_events_to_file(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3], out_path: &str ) {
    let out_img = self.render_events(events, rising_pix, falling_pix);
    out_img.save(out_path).expect("Couldn't save");
  }

  /// Render events as corners to an image result
  pub fn render_corners(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3] ) -> RgbImage {
    let mut out_img =   RgbImage::new(self.n_pixel_cols , self.n_pixel_rows);

    for evt in events {
      let px = match evt.polarity {
        1 => image::Rgb(*rising_pix),
        0 | _ => image::Rgb(*falling_pix),
      };

      // place a cross at every corner
      drawing::draw_cross_mut(&mut out_img, px, evt.col as i32, evt.row as i32);
    }

    out_img
  }

  /// render events as corners to a file
  pub fn render_corners_to_file(&self, events: &Vec<SaeEvent>, rising_pix:&[u8;3], falling_pix:&[u8;3], out_path: &str ) {
    let out_img = self.render_corners(events, rising_pix, falling_pix);
    out_img.save(out_path).expect("Couldn't save");

  }

}



#[cfg(test)]
mod tests {
  use super::*;
  use std::fs::create_dir_all;
  use std::path::Path;

  pub fn insert_synth_feature(tracker: &mut Box<FeatureTracker>,
                              row: i32, col: i32,
                              timestamp: &mut SaeTime, polarity: u8,
                              desc_val: f32 )  -> Option<SaeEvent> {
    let ndesc:NormDescriptor = [desc_val; NORM_DESCRIPTOR_LEN];
    *timestamp += SYNTH_TIME_INCR;

    let evt = SaeEvent {
      row: row as u16,
      col: col as u16,
      polarity,
      timestamp: *timestamp,
      norm_descriptor: Some(Box::new(ndesc)),
    };

    let parent_opt = tracker.track_feature(&evt);
    parent_opt
  }


  const TIMESCALE: f64 = 1E-6; // 1 microsecond per SaeTick
  const SYNTH_TIME_INCR:SaeTime = 10; // ticks between synthetic features

  /// Generate a series of synthetic events for stimulating the tracker
  pub fn process_worm_race_features(img_w: u32, img_h: u32, render_out: bool) {
    let time_window = (0.1 / TIMESCALE) as SaeTime; //0.1 second
    let ref_time_filter = 0; //(50E-3 / TIMESCALE) as SaeTime; //50ms

    let mut tracker =
      Box::new(FeatureTracker::new(img_w, img_h, time_window, ref_time_filter));
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

    let max_time_delta: SaeTime = 20 * SYNTH_TIME_INCR ;

    for _i in 0..total_frames {
      chunk_count += 1;

      insert_synth_feature(&mut tracker,ctr_row - row_gap, ctr_col, &mut timestamp, 0, 0.2);
      insert_synth_feature(&mut tracker, ctr_row, ctr_col,  &mut timestamp, 0, 0.4);
      insert_synth_feature(&mut tracker,ctr_row + row_gap, ctr_col, &mut timestamp, 0, 0.6);

      insert_synth_feature(&mut tracker,ctr_row - row_gap, ctr_col - col_gap,  &mut timestamp, 1, 0.2);
      insert_synth_feature(&mut tracker,ctr_row, ctr_col - col_gap, &mut timestamp, 1, 0.4);
      insert_synth_feature(&mut tracker,ctr_row + row_gap, ctr_col - col_gap, &mut timestamp, 1, 0.6);

      let horizon = timestamp.max( max_time_delta) - max_time_delta;
      let tracks = tracker.collect_tracks(2, horizon);
      println!("chunk {} tracks_len: {}", chunk_count, tracks.len());

      if chunk_count > 2 {
        assert_eq!(tracks.len(), 6);
      }

      if render_out {
        let out_path= format!("./out/sae_tracks_{:04}.png", chunk_count);
        tracker.render_tracks_to_file(horizon, &out_path);
      }

      ctr_col -= COL_DECR;
    }
  }

  #[test]
  fn test_synthetic_worm_race() {
    let img_w = 160;
    let img_h = 160;
    create_dir_all(Path::new("./out/")).expect("Couldn't create output dir");
    process_worm_race_features(img_w, img_h, true);
  }
}