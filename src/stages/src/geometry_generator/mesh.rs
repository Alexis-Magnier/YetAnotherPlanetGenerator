use std::num::TryFromIntError;

/// Sentinel value used to represent an empty or non-existent halfedge index.
pub const EMPTY: u32 = u32::MAX;

/// A dual mesh graph structure combining Delaunay triangulation and Voronoi diagrams.
///
/// # Element Indexing
/// - `0 <= r < num_regions`
/// - `0 <= s < num_sides`
/// - `0 <= t < num_triangles`
///
/// # Ghost & Boundary Elements
/// Outer boundary triangles and regions are wrapped around a synthetic "ghost region"
/// at `num_regions - 1` to close the topological manifold without special null-checks.
#[derive(Debug, Clone, Default)]
pub struct Mesh {
    /// Number of regions marked on the boundary.
    pub num_boundary_regions: u32,
    /// Number of non-ghost (solid) directed sides in the mesh.
    pub num_solid_sides: u32,
    /// Total number of directed sides, including ghost sides.
    pub num_sides: u32,
    /// Total number of regions, including the ghost region.
    pub num_regions: u32,
    /// Total number of non-ghost (solid) regions.
    pub num_solid_regions: u32,
    /// Total number of triangles, including ghost triangles.
    pub num_triangles: u32,
    /// Total number of non-ghost (solid) triangles.
    pub num_solid_triangles: u32,

    /// Map from directed side `s` to its starting region ID `r`.
    pub triangles: Vec<u32>,
    /// Map from directed side `s` to its opposite halfedge side ID `s`.
    pub halfedges: Vec<u32>,

    /// Index mapping each region `r` to an incoming directed side `s`.
    pub r_in_s: Vec<u32>,
}

impl Mesh {
    pub const EMPTY: u32 = EMPTY;

    /// Converts a side ID `s` to its parent triangle ID `t`.
    ///
    /// # Example
    /// ```rust
    /// # use dual_mesh::Mesh;
    /// assert_eq!(Mesh::s_to_t(7), 2);
    /// ```
    #[inline]
    pub fn s_to_t(s: u32) -> u32 {
        s / 3
    }

    /// Returns the previous directed side ID `s` within the same triangle.
    #[inline]
    pub fn s_prev_s(s: u32) -> u32 {
        if s % 3 == 0 { s + 2 } else { s - 1 }
    }

    /// Returns the next directed side ID `s` within the same triangle.
    #[inline]
    pub fn s_next_s(s: u32) -> u32 {
        if s % 3 == 2 { s - 2 } else { s + 1 }
    }

    /// Creates a new `Mesh` from partial mesh data and computes derived attributes.
    ///
    /// # Parameters
    /// - `num_boundary_regions`: Number of boundary regions present in the layout.
    /// - `num_solid_sides`: Number of non-ghost sides.
    /// - `num_regions`: Number of regions
    /// - `triangles`: Halfedge triangle corner map.
    /// - `halfedges`: Halfedge twin map.
    pub fn new(
        num_boundary_regions: u32,
        num_solid_sides: u32,
        num_regions: u32,
        triangles: Vec<u32>,
        halfedges: Vec<u32>,
    ) -> Self {
        let mut mesh = Self {
            num_boundary_regions,
            num_solid_sides,
            num_sides: 0,
            num_regions,
            num_solid_regions: 0,
            num_triangles: 0,
            num_solid_triangles: 0,
            triangles,
            halfedges,
            r_in_s: Vec::new(),
        };
        mesh.update_internal();
        mesh
    }

    /// Constructs a `Mesh` directly from raw Delaunator output without boundary regions.
    ///
    /// # Parameters
    /// - `points`: Flat vector of region positions `[x, y]`.
    /// - `triangles`: Triangulation array produced by Delaunator.
    /// - `halfedges`: Halfedge array produced by Delaunator.
    pub fn from_delaunator(
        num_regions: u32,
        triangles: Vec<u32>,
        halfedges: Vec<u32>,
    ) -> Result<Self, String> {

        let num_solid_sides: u32 = triangles.len()
            .try_into()
            .map_err(|err: TryFromIntError| err.to_string())?;

        Ok(Self::new(0, num_solid_sides, num_regions, triangles, halfedges))
    }

    /// Updates existing mesh arrays with new Delaunator data and recalculates internal state.
    pub fn update(&mut self, triangles: Vec<u32>, halfedges: Vec<u32>) {
        self.triangles = triangles;
        self.halfedges = halfedges;
        self.update_internal();
    }

    /// Recalculates element counts, region-to-side mappings, and triangle centroid positions.
    fn update_internal(&mut self) {
        self.num_sides = self.triangles.len() as u32;
        self.num_solid_regions = self.num_regions.saturating_sub(1);
        self.num_triangles = self.num_sides / 3;
        self.num_solid_triangles = self.num_solid_sides / 3;

        // Construct index mapping region ID -> incoming side ID
        self.r_in_s.clear();
        self.r_in_s.resize(self.num_regions as usize, EMPTY);
        
        for s in 0..self.triangles.len() as u32 {
            let endpoint = self.triangles[Self::s_next_s(s) as usize];
            if self.r_in_s[endpoint as usize] == EMPTY || self.halfedges[s as usize] == EMPTY {
                self.r_in_s[endpoint as usize] = s;
            }
        }
    }

    // --- Accessors ---
    
    /// Returns the starting region ID `r` for side `s`.
    #[inline] pub fn s_begin_r(&self, s: u32) -> u32 { self.triangles[s as usize] }

    /// Returns the ending region ID `r` for side `s`.
    #[inline] pub fn s_end_r(&self, s: u32) -> u32 { self.triangles[Self::s_next_s(s) as usize] }

    /// Returns the interior triangle ID `t` containing side `s`.
    #[inline] pub fn s_inner_t(&self, s: u32) -> u32 { Self::s_to_t(s) }

    /// Returns the exterior triangle ID `t` adjacent to side `s` across its halfedge twin.
    #[inline] pub fn s_outer_t(&self, s: u32) -> u32 { Self::s_to_t(self.halfedges[s as usize]) }

    /// Returns the twin/opposite side ID `s` for side `s`.
    #[inline] pub fn s_opposite_s(&self, s: u32) -> u32 { self.halfedges[s as usize] }

    // --- Triangle Circulators ---

    /// Returns the 3 directed side IDs belonging to triangle `t`.
    #[inline] pub fn t_circulate_s(&self, t: u32) -> [u32; 3] { [3 * t, 3 * t + 1, 3 * t + 2] }

    /// Returns the 3 region IDs forming the corners of triangle `t`.
    #[inline] pub fn t_circulate_r(&self, t: u32) -> [u32; 3] {
        [
            self.triangles[(3 * t) as usize],
            self.triangles[(3 * t + 1) as usize],
            self.triangles[(3 * t + 2) as usize]]
    }

    /// Returns the 3 neighboring triangle IDs adjacent to triangle `t`.
    #[inline] pub fn t_circulate_t(&self, t: u32) -> [u32; 3] {
        [self.s_outer_t(3 * t), self.s_outer_t(3 * t + 1), self.s_outer_t(3 * t + 2)]
    }

    // --- Region Circulators (Zero Allocation) ---

    /// Returns a zero-allocation iterator over all directed side IDs `s` surrounding region `r`.
    #[inline]
    pub fn r_circulate_s(&self, r: u32) -> RegionCirculator<'_, fn(&Mesh, u32) -> u32> {
        self.create_circulator(r, |mesh, incoming| mesh.halfedges[incoming as usize])
    }

    /// Returns a zero-allocation iterator over all neighboring region IDs `r` surrounding region `r`.
    #[inline]
    pub fn r_circulate_r(&self, r: u32) -> RegionCirculator<'_, fn(&Mesh, u32) -> u32> {
        self.create_circulator(r, |mesh, incoming| mesh.s_begin_r(incoming))
    }

    /// Returns a zero-allocation iterator over all neighboring triangle IDs `t` surrounding region `r`.
    #[inline]
    pub fn r_circulate_t(&self, r: u32) -> RegionCirculator<'_, fn(&Mesh, u32) -> u32> {
        self.create_circulator(r, |_mesh, incoming| Mesh::s_to_t(incoming))
    }

    #[inline]
    fn create_circulator<F>(&self, r: u32, extract: F) -> RegionCirculator<'_, F>
    where
        F: Fn(&Mesh, u32) -> u32,
    {
        let start_s = self.r_in_s[r as usize];
        RegionCirculator {
            mesh: self,
            start_s,
            incoming: start_s,
            done: start_s == EMPTY,
            extract,
        }
    }

    // --- Ghost & Boundary Predicates ---

    /// Returns the ID of the synthetic ghost region.
    #[inline] pub fn ghost_r(&self) -> u32 { self.num_regions - 1 }

    /// Checks if side `s` is a ghost side.
    #[inline] pub fn s_ghost(&self, s: u32) -> bool { s >= self.num_solid_sides }

    /// Checks if region `r` is the ghost region.
    #[inline] pub fn r_ghost(&self, r: u32) -> bool { r == self.num_regions - 1 }

    /// Checks if triangle `t` is a ghost triangle.
    #[inline] pub fn t_ghost(&self, t: u32) -> bool { self.s_ghost(3 * t) }

    /// Checks if side `s` lies on the outer boundary.
    #[inline] pub fn s_boundary(&self, s: u32) -> bool { self.s_ghost(s) && (s % 3 == 0) }

    /// Checks if region `r` is a designated boundary region.
    #[inline] pub fn r_boundary(&self, r: u32) -> bool { r < self.num_boundary_regions }
}


/// Constructed via [`Mesh::r_circulate_s`], [`Mesh::r_circulate_r`], or [`Mesh::r_circulate_t`].
pub struct RegionCirculator<'a, F>
where
    F: Fn(&Mesh, u32) -> u32,
{
    mesh: &'a Mesh,
    start_s: u32,
    incoming: u32,
    done: bool,
    extract: F,
}

impl<'a, F> Iterator for RegionCirculator<'a, F>
where
    F: Fn(&Mesh, u32) -> u32,
{
    type Item = u32;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.done || self.incoming == EMPTY {
            return None;
        }

        let val = (self.extract)(self.mesh, self.incoming);

        let outgoing = Mesh::s_next_s(self.incoming);
        self.incoming = self.mesh.halfedges[outgoing as usize];

        if self.incoming == self.start_s {
            self.done = true;
        }

        Some(val)
    }
}