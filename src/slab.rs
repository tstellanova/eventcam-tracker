
use arcstar::sae_types::*;
use arrayvec::ArrayVec;
use rand::{Rng};


//TODO move some of these constants to a configuration file?

/// Maximum rectilinear taxicab/manhattan spatial distance between center point and matching neighbor
/// (see related "dconn" from arcstar tracking paper, reduced to 5x5 patch in ACE paper)
pub const MAX_RL_SPATIAL_DISTANCE: u32 = 4;

const PATCH_SQUARE_DIM: usize = 5; //5x5 grid around center point
const MAX_PATCH_NEIGHBORS: usize = (PATCH_SQUARE_DIM*PATCH_SQUARE_DIM) - 1;

/// max number of pixels in each 2D dimension
const MAX_PIXEL_DIM: usize = 1024;

/// maximum "tree depth", similar to "ρmax"/"ρthresh" from ACE tracker paper
const TREE_DEPTH: usize = 5;


/// a slab is Surface of Active Events with historical depth
type FeatureChain =  ArrayVec<[SaeEvent; TREE_DEPTH]>;
type SlabRow = ArrayVec<[SlabCell; MAX_PIXEL_DIM]>;
type SlabMatrix = Vec<SlabRow>;

/// row, column for cell in slab
type SlabCellId = (usize, usize);


#[derive( Clone, Debug, PartialEq)]
pub struct FeatureTrack {
  /// stores a chain of features in the track
  pub chain: FeatureChain,
  /// a color used for visually representing this track
  pub color: [u8; 3],
}

type SlabCell = FeatureTrack;

impl std::default::Default for FeatureTrack {
  fn default() -> Self {
    FeatureTrack {
      chain: Default::default(),
      color: [0; 3],
    }
  }
}

impl FeatureTrack {
  pub fn new() -> Self {
    Self::default()
  }
}


#[derive( Clone, Debug, PartialEq)]
pub struct FeatureLocator {
  cell_id: SlabCellId,
  pub event: SaeEvent,
}

impl std::default::Default for FeatureLocator {
  fn default() -> Self {
    FeatureLocator {
      cell_id: (0, 0),
      event: Default::default()
    }
  }
}

pub struct SlabStore {
  mx: SlabMatrix,
}


impl SlabStore {
  pub fn new() -> SlabStore {

    let mut le_row = SlabRow::new();
    for _i in 0..MAX_PIXEL_DIM {
      le_row.push(SlabCell::default());
    }

    let mut le_matrix = SlabMatrix::with_capacity(MAX_PIXEL_DIM);
    for _i in 0..MAX_PIXEL_DIM {
      le_matrix.push(le_row.clone());
    }

    let inst = SlabStore {
      mx: le_matrix,
    };

    inst
  }

  /// Get a track of features for a given location, iff they're fresher than the time horizon
  pub fn track_for_point(&self, row: u16, col: u16, time_horizon: SaeTime) -> Option<FeatureTrack> {

    let slab_cell: &SlabCell = &self.mx[row as usize][col as usize];

    if slab_cell.chain.len() > 1 &&
      slab_cell.chain[0].timestamp > time_horizon   {
      let mut track: FeatureTrack = FeatureTrack::new();
      track.color = slab_cell.color;

      for feature in slab_cell.chain.iter() {
        if feature.timestamp >= time_horizon {
          track.chain.push(feature.clone());
        }
        else { break;} //no more fresh features
      }

      return Some(track)
    }

    None
  }


  //figure 5 from the paper is the closest thing that comes to my mind.
  // Keep in mind that the core structure is really similar to
  // the tracker described in the Arc* paper
  // however, new rules for matching (with a pseudo-descriptor) and tree-growth are proposed.

  /// find best parent node based on spatial and temporal distance
  fn find_closest_and_freshest_neighbor(
    evt: &SaeEvent,
    neighbors: &Vec<FeatureLocator>
  ) -> Option<FeatureLocator> {
    let mut best_neighbor: FeatureLocator = FeatureLocator::default();
    let mut best_neighbor_likeness = 0.0;
    let mut best_neighbor_time: SaeTime = 0;
    let mut best_neighbor_dist= u32::max_value();


    for nb in neighbors {
      let prior_evt = &nb.event;
      if prior_evt.polarity == evt.polarity {
        let spatial_dist = evt.spatial_rl_dist(&prior_evt);
        if spatial_dist <= MAX_RL_SPATIAL_DISTANCE {
          let likeness = evt.likeness(&prior_evt);

          if likeness > 0.5f32 { //dmax = 0.5 from ACE tracker paper
            //prefer the most-like neighbors
            if likeness > best_neighbor_likeness {
              best_neighbor = nb.clone();
              best_neighbor_likeness = likeness;
              best_neighbor_dist = spatial_dist;
              best_neighbor_time = prior_evt.timestamp;
            }
            else if likeness == best_neighbor_likeness {
              // if likness is the same, prefer spatial, then temporal closer neighbors
              if spatial_dist < best_neighbor_dist {
                best_neighbor = nb.clone();
                best_neighbor_dist = spatial_dist;
                best_neighbor_time = prior_evt.timestamp;
              }
              else if spatial_dist == best_neighbor_dist {
                if prior_evt.timestamp > best_neighbor_time {
                  best_neighbor = nb.clone();
                  best_neighbor_dist = spatial_dist;
                  best_neighbor_time = prior_evt.timestamp;
                }
              }
            }

          }
        }
      }
    }

    if best_neighbor_likeness > 0.0 {
      Some(best_neighbor)
    }
    else {
      None
    }

  }

  /// Return the best matching parent cell or none
  fn find_best_parent_cell(&self, evt: &SaeEvent, time_horizon: SaeTime) -> Option<FeatureLocator> {
    let (leaf_neighbors, nonleaf_neighbors) = self.collect_neighbors(evt, time_horizon);
//    println!("leaf_neighbors {} nonleaf_neighbors {}", leaf_neighbors.len(), nonleaf_neighbors.len());

    let mut best_parent_opt:Option<FeatureLocator> =
      Self::find_closest_and_freshest_neighbor(evt, &leaf_neighbors);

    if best_parent_opt.is_none() {
      //search the nonleaf nodes
      best_parent_opt = Self::find_closest_and_freshest_neighbor(evt, &nonleaf_neighbors);
    }

    best_parent_opt
  }

  /// Collect relevant neighbors in the slab for the given event
  /// return leaf, non-leaf neighbor lists
  fn collect_neighbors(&self, evt: &SaeEvent, time_horizon: SaeTime) -> (Vec<FeatureLocator>, Vec<FeatureLocator>) {
    let irow = evt.row as i32;
    let icol = evt.col as i32;
    let mut leaf_neighbors:Vec<FeatureLocator> = Vec::with_capacity(MAX_PATCH_NEIGHBORS);
    let mut nonleaf_neighbors:Vec<FeatureLocator> = Vec::with_capacity(MAX_PATCH_NEIGHBORS);

    let mut x: i32 = 1; //skip 0,0 since that's not a _neighbor_ for matching purposes?
    let mut y: i32 = 0;

    const MAX_IDX:usize = (PATCH_SQUARE_DIM*PATCH_SQUARE_DIM -1 );//one accounts for (0,0) skip
    for _i in 0..MAX_IDX {
      let row: usize = (irow + y) as usize;
      let col: usize = (icol + x) as usize;
      let slab_cell_id: SlabCellId = (row, col);

      //println!("{} ( {}, {} ) ... [{}, {}]", _i, x, y, row, col);
      let cell = &self.mx[row][col];
      let cell_len = cell.chain.len();
      if cell_len > 0 {
        if cell.chain.len() == 1 {
          let node = &cell.chain[0];
          if node.timestamp >= time_horizon &&
            node.timestamp < evt.timestamp  //assume that evt is the freshest
          {
            leaf_neighbors.push(FeatureLocator {
              cell_id: slab_cell_id,
              event: node.clone()
            });
          }
        }
        else {
          //roll up all the kids in the cell
          for node in cell.chain.iter() {
            if node.timestamp >= time_horizon {
              nonleaf_neighbors.push(FeatureLocator {
                cell_id: slab_cell_id,
                event: node.clone()
              });
            } else {
              //TODO else mark for deletion??
              //println!("time {} < {}",node.timestamp, time_horizon);
            }
          }
        }
      }

      if (x.abs() <= y.abs()) && ((x >= 0) || (x != y)) {
        x += if y >= 0 { 1 } else { -1 };
      }
      else {
        y += if x >= 0 { -1 } else { 1 };
      }
    }

    (leaf_neighbors, nonleaf_neighbors)
  }


  /// insert a feature in the slab
  fn push_feature_to_slab(&mut self, evt: &SaeEvent, parent_opt: Option<FeatureLocator>, time_horizon: SaeTime) {
    let mut out_cell:SlabCell = SlabCell::new();
    out_cell.chain.push(evt.clone());
    out_cell.color = rand::thread_rng().gen::<[u8;3]>();

    if parent_opt.is_some() {
      let parent_loc = parent_opt.unwrap();
      let match_cell_id = parent_loc.cell_id;
      let parent_evt = &parent_loc.event;
      //copy over the chain from the matched cell,
      //starting from the actual parent event!
      let in_cell:&mut SlabCell = &mut self.mx[match_cell_id.0][match_cell_id.1];
      let parent_is_leaf = in_cell.chain.len() == 1;

      //since we've already inserted the freshest event, need to ensure
      //we don't overflow the cell
      let max_idx = in_cell.chain.len().min(TREE_DEPTH -1);

      let mut copy_started = false;
      //copy over the older chain of events from the parent chain
      for i in 0..max_idx  {
        let prev_evt = &in_cell.chain[i];
        if !copy_started {
          //wait until we see the actual parent event before starting to copy
          if (prev_evt.timestamp == parent_evt.timestamp) &&
            (prev_evt.row == parent_evt.row) &&
            (prev_evt.col == parent_evt.col) {
            copy_started = true;
          }
        }

        if copy_started {
          if prev_evt.timestamp >= time_horizon {
            out_cell.chain.push(prev_evt.clone());
          }
          else {
            //the chain contains a stale event -- stop following the chain
            break;
          }
        }
      }

      // if the matched cell is a leaf, also push the new event there
      //TODO add unit test for this odd double-insert of the new event in the slab
      if parent_is_leaf {
        in_cell.chain.insert(0, evt.clone());
      }
      else {
        out_cell.color = in_cell.color;
      }
    }

    self.mx[evt.row as usize][evt.col as usize] = out_cell;

  }

  /// returns any matching parent event
  pub fn add_and_match_feature(&mut self, evt: &SaeEvent, time_horizon: SaeTime) -> Option<SaeEvent> {
    let parent_opt = self.find_best_parent_cell(evt, time_horizon);

    self.push_feature_to_slab(evt, parent_opt.clone(), time_horizon);
    if parent_opt.is_some() {
      Some(parent_opt.unwrap().event)
    }
    else {
      None
    }
  }

}

#[cfg(test)]
mod tests {
  use super::*;

  fn generate_test_event() -> SaeEvent {
    generate_test_event_with_desc(&[0.5f32; NORM_DESCRIPTOR_LEN])
  }

  fn generate_test_event_with_desc(desc: &[f32; NORM_DESCRIPTOR_LEN]) -> SaeEvent {
    let mut evt = SaeEvent::new();
    evt.norm_descriptor = Some(Box::new(*desc));
    evt
  }

  #[test]
  fn test_spiral_insert() {
    // Test whether we can insert features in a spiral starting from a center point,
    // and look to see whether the match tracking is accurate.
    // Note that all of the features inserted in this test have the same descriptor
    let mut store = SlabStore::new();
    let time_horizon: SaeTime = 0;

    let irow: i32 = 320;
    let icol: i32 = 320;
    const SQUARE_DIM: usize = 9;

    let mut prior_node: Option<SaeEvent> = None;

    let mut first_node = SaeEvent::new();
    first_node.timestamp = 0 as SaeTime;
    first_node.row = irow as u16;
    first_node.col = icol as u16;
    let check = store.add_and_match_feature(&first_node, time_horizon);
    assert_eq!(true, check.is_none()); //should be no matches yet

    const MAX_IDX: usize = SQUARE_DIM * SQUARE_DIM;
    let mut x: i32 = 0; // insert at 0,0 first, but don't try to match
    let mut y: i32 = 0;
    for i in 0..MAX_IDX {
      let row: usize = (irow + y) as usize;
      let col: usize = (icol + x) as usize;

      let mut node:SaeEvent = generate_test_event();
      node.timestamp = (i * 100)  as SaeTime;
      node.row = row as u16;
      node.col = col as u16;
      node.norm_descriptor = Some(Box::new([0.5f32; NORM_DESCRIPTOR_LEN]) );
      println!("insert {} ( {}, {} ) -> [{}, {}, {}]", i, x, y, row, col, node.timestamp);

      let check = store.add_and_match_feature(&node, time_horizon);
      if prior_node.is_some() {
        assert_eq!(true, check.is_some());
        //because we're inserting in a spiral, the best match should be prior insertion
        let prior_evt = prior_node.unwrap();
        let check_evt = check.unwrap();
        println!("best for {} check {} s/b {} ",node.timestamp, check_evt.timestamp, prior_evt.timestamp);
        assert_eq!(prior_evt.timestamp, check_evt.timestamp);
      }
      prior_node = Some(node);

      if (x.abs() <= y.abs()) && ((x >= 0) || (x != y)) {
        x += if y >= 0 { 1 } else { -1 };
      } else {
        y += if x >= 0 { -1 } else { 1 };
      }
    }

    assert_eq!(true, prior_node.is_some());
    let last_evt = prior_node.unwrap();
    let track_opt = store.track_for_point(last_evt.row, last_evt.col, time_horizon);
    assert_eq!(track_opt.is_some(), true);
    assert_eq!(track_opt.unwrap().chain.len(), TREE_DEPTH);
  }

  #[test]
  fn test_distant_feature_matching() {
    let mut store = SlabStore::new();
    let mut node1 = generate_test_event();
    let time_horizon: SaeTime = 0;
    node1.timestamp = 100;
    node1.row = 320;
    node1.col = 320;
    let check1 = store.add_and_match_feature(&node1, time_horizon);
    assert_eq!(true, check1.is_none());

    //insert a event that is too far away to be considered a match for the previous event
    let mut node2 = generate_test_event();
    node2.timestamp = 200;
    node2.row = node1.row + 10*MAX_RL_SPATIAL_DISTANCE as u16;
    node2.col = node1.col + 10*MAX_RL_SPATIAL_DISTANCE as u16;
    let check2 = store.add_and_match_feature(&node2, time_horizon);
    assert_eq!(true, check2.is_none());

    //insert event that is a match to the previous event
    let mut node3 = generate_test_event();
    node3.timestamp = 300;
    node3.row = node2.row + (MAX_RL_SPATIAL_DISTANCE/2) as u16;
    node3.col = node2.col + (MAX_RL_SPATIAL_DISTANCE/2) as u16;
    let check3 = store.add_and_match_feature(&node3, time_horizon);
    assert_eq!(true, check3.is_some());
    assert_eq!(check3.unwrap().timestamp, node2.timestamp);

    //verify chain is correct
    let track_opt = store.track_for_point(node3.row, node3.col, time_horizon);
    assert_eq!(track_opt.is_some(), true);
    let chain = track_opt.unwrap().chain;
    assert_eq!(true,  chain.len() == 2);
    assert_eq!(chain[0].timestamp, node3.timestamp);
    assert_eq!(chain[1].timestamp, node2.timestamp);

  }

  #[test]
  fn test_spatial_feature_matching() {
    // note that all of the features inserted in this test have the same descriptor
    let mut store = SlabStore::new();
    let mut node1 = generate_test_event();
    let mut time_horizon: SaeTime = 0;
    node1.timestamp = 100;
    node1.row = 320;
    node1.col = 320;
    store.push_feature_to_slab(&node1, None, time_horizon);
    let track_opt = store.track_for_point(node1.row, node1.col, time_horizon);
    assert_eq!(track_opt.is_none(), true);//no "track", just a single feature

    let parent_loc:FeatureLocator = FeatureLocator {
      cell_id: (node1.row as usize, node1.col as usize),
      event: node1.clone(),
    };

    let mut node2 = node1.clone();
    node2.timestamp = 200;
    store.push_feature_to_slab(&node2, Some(parent_loc), time_horizon);
    //verify that node2 is part of a chain of events including node1
    let track_opt = store.track_for_point(node2.row, node2.col, time_horizon);
    assert_eq!(track_opt.unwrap().chain.len(), 2);

    //add a new spatially close feature that should attach to the existing chain
    let mut node3 = node1.clone();
    node3.col = 322; //move slightly, max distance limited by PATCH_SQUARE_DIM
    node3.timestamp = 300;
    let parent_opt = store.add_and_match_feature(&node3, time_horizon);
    assert_eq!(parent_opt.is_some(), true);
    assert_eq!(parent_opt.unwrap().timestamp, node2.timestamp);

    let track_opt = store.track_for_point(node3.row, node3.col, time_horizon);
    assert_eq!(track_opt.unwrap().chain.len(), 3);

    //add another feature that is nearby but is much newer
    let mut node4 = node1.clone();
    node4.col = 324;  //move slightly, max distance limited by PATCH_SQUARE_DIM
    node4.timestamp = 600;
    time_horizon = 200; //this should cause node1 to be ignored in chain result
    let parent_opt = store.add_and_match_feature(&node4, time_horizon);
    assert_eq!(parent_opt.is_some(), true);
    let track_opt = store.track_for_point(node4.row, node4.col, time_horizon);
    assert_eq!(track_opt.is_some(), true);
    let chain = track_opt.unwrap().chain;
    assert_eq!(chain.len(), 3); // 4, 3, 2
    assert_eq!(chain[0].timestamp, node4.timestamp);

  }

  #[test]
  fn test_likeness_feature_matching() {
    let mut store = SlabStore::new();
    let mut node1 = generate_test_event();
    let time_horizon: SaeTime = 0;
    node1.timestamp = 100;
    node1.row = 320;
    node1.col = 320;

    store.push_feature_to_slab(&node1, None, time_horizon);
    let track_opt = store.track_for_point(node1.row, node1.col, time_horizon);
    assert_eq!(track_opt.is_none(), true); //there's no "chain" just a single feature

    let mut node2 = generate_test_event_with_desc(&[0.7f32; NORM_DESCRIPTOR_LEN]);
    node2.timestamp = 200;
    node2.row = node1.row;
    node2.col = 322; //move slightly
    let parent_opt = store.add_and_match_feature(&node2, time_horizon);
    assert_eq!(parent_opt.is_some(), true);
    assert_eq!(parent_opt.unwrap().timestamp, node1.timestamp);

    let mut node3 = node1.clone();
    node3.col = 321; //move slightly, max distance limited by PATCH_SQUARE_DIM
    node3.timestamp = 300;

    //this node is same distance from both node1 and node2, but should likeness-match with node1
    let parent_opt = store.add_and_match_feature(&node3, time_horizon);
    assert_eq!(parent_opt.is_some(), true);
    assert_eq!(parent_opt.unwrap().timestamp, node1.timestamp);

  }

}

