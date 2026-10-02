use core;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use glam::Vec3;

use super::{Mesh, GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct CellGenerator {
    pub seed: u64,
    pub cell_count: u32,
    pub jitter: f32,
}


struct TriangulationResult{
    pub triangles: Vec<u32>,
    pub halfedges: Vec<u32>
}

impl CellGenerator{
    fn make_rng(&self) -> impl Rng {
        StdRng::seed_from_u64(self.seed)
    }

    fn generate_fibbonacy_points(&self, mut rng: impl Rng) -> Vec<Vec3>{
        let n = self.cell_count as usize;
        let jitter = self.jitter as f64;

        if n == 1 {
            return vec![Vec3::new(0.0, 1.0, 0.0)];
        }

        let mut points: Vec<Vec3> = Vec::with_capacity(n);

        const PI:f64 = core::f64::consts::PI;

        let s = 3.6 / (n as f64).sqrt();
        let dlong = PI * (3.- f64::sqrt(5.));  /* ~2.39996323 */
        let dz = 2.0 / (n as f64);

        let mut long:f64 = 0.0;
        let mut z = 1. - dz/2.;

        for _ in 0..n {

            let r = (1. - z*z).sqrt();
            let mut lat = z.asin();
            let mut lon = long;

            let rand_lat = rng.next_u32() as f64 / u32::MAX as f64;
            let rand_long = rng.next_u32() as f64 / u32::MAX as f64;

            lat += jitter * rand_lat * (lat - f64::asin(f64::max(-1., z - dz * 2. * PI * r / s)));
            lon += jitter * rand_long * s/r;

            points.push(Vec3{
                x: (lat.cos() * lon.cos()) as f32,
                y: (lat.cos() * lon.sin()) as f32,
                z: (lat.sin()) as f32,
            });
        
            z -= dz;
            long += dlong;
        }
        points
    }

    fn triangulate(&self, points: &mut Vec<Vec3>) -> TriangulationResult{
        use delaunator::{EMPTY, Point, Triangulation, triangulate};

        fn stereographic_projection(points: &Vec<Vec3>) -> Vec<Point> {
            points.iter().map(|point|{
                Point{
                    x: (point.x / (1. - point.z)) as f64,
                    y: (point.y / (1. - point.z)) as f64
                }
            }).collect()
        }

        fn add_pole(pole_id: usize, triangulation: &mut Triangulation){

            let n = triangulation.triangles.len();

            fn next(i: usize) -> usize{
                if (i % 3) == 2 {
                    i - 2
                } else {
                    i + 1
                }
            }

            let mut num_unpaired: usize = 0;
            let mut first_unpaired: usize = 0;

            let mut point_to_side = vec![0; n];

            for s in 0..n {
                if triangulation.halfedges[s] == EMPTY {
                    num_unpaired += 1;
                    point_to_side[triangulation.triangles[s]] = s;
                    first_unpaired = s;
                }
            }

            let total_capacity = n + 3 * num_unpaired;

            triangulation.triangles.resize(total_capacity, 0);
            triangulation.halfedges.resize(total_capacity, 0);

            let mut s = first_unpaired;
            for i in 0..num_unpaired {

                let ns = n + 3 * i;

                triangulation.halfedges[s] = ns;
                triangulation.halfedges[ns] = s;
                triangulation.triangles[ns]     = triangulation.triangles[next(s)];
                triangulation.triangles[ns + 1] = triangulation.triangles[s];
                triangulation.triangles[ns + 2] = pole_id;

                let k = n + (3 * i + 4) % (3 * num_unpaired);

                triangulation.halfedges[ns + 2]  = k;
                triangulation.halfedges[k]       = ns + 2;

                s = point_to_side[triangulation.triangles[next(s)]]
            }
        }

        let flat = stereographic_projection(&points);
        let mut result = triangulate(&flat);

        add_pole(points.len(), &mut result);
        
        points.push(Vec3::new(0.0, 0.0, 1.0));
        
        TriangulationResult{
            triangles: result.triangles.iter().map(|t| *t as u32 ).collect(),
            halfedges: result.halfedges.iter().map(|h| *h as u32 ).collect(), 
        }
    }
    
}

impl GeneratorStage for CellGenerator{
    fn name(&self) -> &'static str {
        "Cell generator"
    }

    fn run(&self, data: &mut PlanetData) -> anyhow::Result<()> {

        let rng = self.make_rng();
        let mut points = self.generate_fibbonacy_points(rng);

        let triangulation = self.triangulate(&mut points);


        let mesh = Mesh::from_delaunator(
            self.cell_count + 1,
            triangulation.triangles,
            triangulation.halfedges
        )
            .map_err(|e| anyhow::anyhow!(e))?;

        data.mesh = mesh;
        data.cells.position = points;
        data.cells.count = self.cell_count+1;

        Ok(())
    }
}