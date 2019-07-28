
use arcstar::sae_types::*;

use arrayvec::ArrayVec;



//TODO move some of these constants to a configuration file?

/// Maximum rectilinear taxicab/manhattan spatial distance between center point and matching neighbor
/// (see related "dconn" from arcstar tracking paper, reduced to 5x5 patch in ACE paper)
const MAX_RL_SPATIAL_DISTANCE: u32 = 4;

const PATCH_SQUARE_DIM: usize = 5; //5x5 grid around center point

/// max number of pixels in each 2D dimension
const MAX_PIXEL_DIM: usize = 1024;

/// maximum "tree depth", similar to "ρmax"/"ρthresh" from ACE tracker paper
pub const SLAB_DEPTH: usize = 5;


/// a slab is Surface of Active Events with historical depth
type SlabCell =  ArrayVec<[SaeEvent;SLAB_DEPTH]>;
type SlabRow = ArrayVec<[SlabCell; MAX_PIXEL_DIM]>;
type SlabMatrix = Vec<SlabRow>;

type SlabCellId = (usize, usize);

#[derive( Clone, Debug, PartialEq)]
pub struct FeatureLocator {
  cell_id: SlabCellId,
  pub event: SaeEvent,
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

  pub fn chain_for_point(&self, row: u16, col: u16, time_horizon: SaeTime) -> Vec<SaeEvent> {
    let mut res:Vec<SaeEvent> = vec!();
    let cell = &self.mx[row as usize][col as usize];
    for item in cell {
      if item.timestamp >= time_horizon {
        res.push(item.clone());
      }
    }

    res
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
    let mut best_parent_opt:Option<FeatureLocator> = None;
    let mut best_parent_time: SaeTime = 0;
    let mut best_parent_distance= u32::max_value();

    for nb in neighbors {
      let prior_evt = &nb.event;
      if prior_evt.polarity == evt.polarity {
        let spatial_dist = evt.spatial_rl_dist(&prior_evt);
        if spatial_dist <= MAX_RL_SPATIAL_DISTANCE {
          let likeness = evt.likeness(&prior_evt);
          if likeness > 0.5f32 { //dmax = 0.5 from ACE tracker paper
            if spatial_dist < best_parent_distance {
              best_parent_distance = spatial_dist;
              best_parent_time = prior_evt.timestamp;
              best_parent_opt = Some(nb.clone());
            }
            else if spatial_dist == best_parent_distance {
              if prior_evt.timestamp > best_parent_time {
                best_parent_distance = spatial_dist;
                best_parent_time = prior_evt.timestamp;
                best_parent_opt = Some(nb.clone());
              }
            }
          }
        }
      }
    }

    best_parent_opt
  }

  /// Return the best matching parent cell or none
  fn find_best_parent_cell(&self, evt: &SaeEvent, time_horizon: SaeTime) -> Option<FeatureLocator> {
    let (leaf_neighbors, nonleaf_neighbors) = self.collect_neighbors(evt, time_horizon);
    //println!("leaf_neighbors {} nonleaf_neighbors {}", leaf_neighbors.len(), nonleaf_neighbors.len());

    let mut best_parent_opt:Option<FeatureLocator> =
      Self::find_closest_and_freshest_neighbor(evt, &leaf_neighbors);

    if best_parent_opt.is_none() {
      //search the nonleaf nodes
      best_parent_opt = Self::find_closest_and_freshest_neighbor(evt, &nonleaf_neighbors);
    }

    best_parent_opt
  }

  /// Collect relevant neighbors in the forest for the given event
  /// return leaf, non-leaf neighbor lists
  fn collect_neighbors(&self, evt: &SaeEvent, time_horizon: SaeTime) -> (Vec<FeatureLocator>, Vec<FeatureLocator>) {
    let irow = evt.row as i32;
    let icol = evt.col as i32;
    let mut leaf_neighbors:Vec<FeatureLocator> = vec!();
    let mut nonleaf_neighbors:Vec<FeatureLocator> = vec!();
//    let time_horizon: SaeTime = evt.timestamp.max(FORGETTING_TIME) - FORGETTING_TIME;

    let mut x: i32 = 1; //skip 0,0 since that's not a _neighbor_ for matching purposes?
    let mut y: i32 = 0;

    const MAX_IDX:usize = (PATCH_SQUARE_DIM*PATCH_SQUARE_DIM -1 );//one accounts for (0,0) skip
    for _i in 0..MAX_IDX {
      let row: usize = (irow + y) as usize;
      let col: usize = (icol + x) as usize;
      let slab_cell_id: SlabCellId = (row, col);

      //println!("{} ( {}, {} ) ... [{}, {}]", i, x, y, row, col);
      let cell = &self.mx[row][col];
      if cell.len() == 1 {
        let node = &cell[0];
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
        for node in cell {
          if node.timestamp >= time_horizon {
            nonleaf_neighbors.push(FeatureLocator {
              cell_id: slab_cell_id,
              event: node.clone()
            });
          }
          else {
            //TODO else mark for deletion??
            //println!("time {} < {}",node.timestamp, time_horizon);
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


  /// returns any matching parent event
  pub fn add_and_match_feature(&mut self, evt: &SaeEvent, time_horizon: SaeTime) -> Option<FeatureLocator> {
//    let time_horizon: SaeTime = evt.timestamp.max(FORGETTING_TIME) - FORGETTING_TIME;

    let mut matched_feature_opt = None;

    let mut out_cell:SlabCell = SlabCell::new();
    //insert freshest event at the top of the cell
    out_cell.push(evt.clone());


    let parent = self.find_best_parent_cell(evt, time_horizon);
    if parent.is_some() {
      matched_feature_opt = parent.clone();

      let match_cell_id = parent.unwrap().cell_id;
      //copy over the chain from the matched cell
      let in_cell:&mut SlabCell = &mut self.mx[match_cell_id.0][match_cell_id.1];
      let parent_is_leaf = in_cell.len() == 1;

      //since we've already inserted the freshest event, need to ensure
      //we don't overflow the cell
      let max_idx = in_cell.len().min(SLAB_DEPTH-1);

      //copy over the older chain of events from the parent chain
      for i in 0..max_idx  {
        if in_cell[i].timestamp >= time_horizon {
          out_cell.push(in_cell[i].clone());
        }
      }

      // if the matched cell is a leaf, also push the new event there
      if parent_is_leaf {
        in_cell.insert(0, evt.clone());
      }
    }
    else {
      //println!("no parent for {:?}", evt);
    }


    self.mx[evt.row as usize][evt.col as usize] = out_cell;

    matched_feature_opt
  }

}

#[cfg(test)]
mod tests {
  use super::*;

  fn generate_test_event() -> SaeEvent {
    let mut evt = SaeEvent::new();
    evt.norm_descriptor = Some(Box::new([0.5; NORM_DESCRIPTOR_LEN]));
    evt
  }

  #[test]
  fn test_spiral_insert() {
    // Test whether we can insert events in a spiral starting from a center point,
    // and look to see whether the match tracking is accurate.
    let mut store = SlabStore::new();
    let time_horizon: SaeTime = 0;

    let irow: i32 = 320;
    let icol: i32 = 320;
    const SQUARE_DIM: usize = 9;

    let mut prior_node: Option<FeatureLocator> = None;

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
        let prior_evt = prior_node.unwrap().event;
        let check_evt = check.unwrap().event;
        println!("best for {} check {} s/b {} ",node.timestamp, check_evt.timestamp, prior_evt.timestamp);
        assert_eq!(prior_evt.timestamp, check_evt.timestamp);
      }
      prior_node = Some(FeatureLocator {
        cell_id: (node.row as usize, node.col as usize),
        event: node
      });

      if (x.abs() <= y.abs()) && ((x >= 0) || (x != y)) {
        x += if y >= 0 { 1 } else { -1 };
      } else {
        y += if x >= 0 { -1 } else { 1 };
      }
    }

    assert_eq!(true, prior_node.is_some());
    let last_evt = prior_node.unwrap().event;
    let chain = store.chain_for_point(last_evt.row, last_evt.col, time_horizon);
    assert_eq!(chain.len(), SLAB_DEPTH);
  }



  #[test]
  fn test_distant_feature_insertion() {
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
    assert_eq!(check3.unwrap().event.timestamp, node2.timestamp);

    //verify chain is correct
    let chain = store.chain_for_point(node3.row, node3.col, time_horizon);
    assert_eq!(true,  chain.len() == 2);
    assert_eq!(chain[0].timestamp, node3.timestamp);
    assert_eq!(chain[1].timestamp, node2.timestamp);

  }

}

