# Simple BVH based ray tracing kernel for Parallella Epiphany.

Status: just beginning. Nothing working yet.

## Design

* Try to code&data all fit into Epiphany on-chip memory(32KB, 8KB x 4 banks, for each Epiphany core)
  * Render small scene(~65,536 triangles)

## TODO

*  [x] Ray - AABB intersection
*  [ ] Ray - Triangle intersection
*  [ ] BVH Traversal

## Performance

* Ray - AABB intersection: 100 clocks
 
## Note

* On-chip data load/store: 8byte/cycle?
* Core to core data transfer: 1,000~ cycles!
* Core to DRRAM data transfer: 10,000~ cycles!
