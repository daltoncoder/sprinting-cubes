use std::collections::HashMap;

use crate::tables::{CORNER_OFFSETS, EDGE_CORNERS, EDGE_TABLE, TRI_TABLE};

pub struct Mesh {
    pub vertices: Vec<[f32; 3]>,        // Unique vertices
    pub triangles: Vec<[usize; 3]>,     // Indices into vertices
    pub normals: Option<Vec<[f32; 3]>>, // Per-vertex normals
}

impl Mesh {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
            normals: None,
        }
    }
    // Convert 3D coordinates to flat index
    pub fn coord_to_index(x: usize, y: usize, z: usize, dims: (usize, usize, usize)) -> usize {
        let (width, height, _depth) = dims;
        z * height * width + y * width + x
    }

    // Convert flat index back to 3D coordinates
    pub fn index_to_coord(index: usize, dims: (usize, usize, usize)) -> (usize, usize, usize) {
        let (width, height, _depth) = dims;
        let z = index / (height * width);
        let remainder = index % (height * width);
        let y = remainder / width;
        let x = remainder % width;
        (x, y, z)
    }
}

/// Performs marching cubes on a 3D scalar field.
///
/// # Arguments
/// * `points` - Flattened 3D array in Z-Y-X order (z varies slowest)
/// * `dimensions` - (width, height, depth) of the 3D grid
/// * `isolevel` - Surface threshold value
///
/// # Coordinate System
/// * Origin at bottom-left-front corner
/// * X-axis: left to right
/// * Y-axis: bottom to top  
/// * Z-axis: front to back
///
/// # Array Layout
/// ```
/// // For a 3x2x2 grid:
/// let points = [
///     // Z=0 slice (front)
///     // Y=0 row (bottom)
///     points[0], points[1], points[2],  // (0,0,0), (1,0,0), (2,0,0)
///     // Y=1 row (top)  
///     points[3], points[4], points[5],  // (0,1,0), (1,1,0), (2,1,0)
///     
///     // Z=1 slice (back)
///     // Y=0 row (bottom)
///     points[6], points[7], points[8],  // (0,0,1), (1,0,1), (2,0,1)
///     // Y=1 row (top)
///     points[9], points[10], points[11], // (0,1,1), (1,1,1), (2,1,1)
/// ];
/// ```
pub fn marching_cubes(points: &[f32], dimensions: (usize, usize, usize), isolevel: f32) -> Mesh {
    let (width, height, depth) = dimensions;
    // Validate input
    if points.len() != width * height * depth {
        panic!("Points array size doesn't match dimensions");
    }

    if width < 2 || height < 2 || depth < 2 {
        return Mesh::new(); // Need at least 2x2x2 to form cubes
    }

    let mut mesh = Mesh::new();

    let mut vertex_cache: HashMap<(usize, usize, usize, usize), usize> = HashMap::new();

    // Standard iteration order:
    // z varies slowest (outermost loop)
    // y varies middle
    // x varies fastest (innermost loop)
    for z in 0..depth - 1 {
        for y in 0..height - 1 {
            for x in 0..width - 1 {
                // Process point at (x, y, z)
                process_cube(
                    points,
                    dimensions,
                    (x, y, z),
                    isolevel,
                    &mut mesh,
                    &mut vertex_cache,
                );
            }
        }
    }
    mesh
}

fn process_cube(
    points: &[f32],
    dimensions: (usize, usize, usize),
    cube_pos: (usize, usize, usize),
    isolevel: f32,
    mesh: &mut Mesh,
    vertex_cache: &mut HashMap<(usize, usize, usize, usize), usize>,
) {
    let (x, y, z) = cube_pos;

    // Get the 8 corner values of this cube
    let cube_values = [
        get_point_value(points, dimensions, (x, y, z)), // 0: bottom-left-front
        get_point_value(points, dimensions, (x + 1, y, z)), // 1: bottom-right-front
        get_point_value(points, dimensions, (x + 1, y + 1, z)), // 2: top-right-front
        get_point_value(points, dimensions, (x, y + 1, z)), // 3: top-left-front
        get_point_value(points, dimensions, (x, y, z + 1)), // 4: bottom-left-back
        get_point_value(points, dimensions, (x + 1, y, z + 1)), // 5: bottom-right-back
        get_point_value(points, dimensions, (x + 1, y + 1, z + 1)), // 6: top-right-back
        get_point_value(points, dimensions, (x, y + 1, z + 1)), // 7: top-left-back
    ];

    // Calculate cube configuration index
    let mut cube_index = 0;
    for i in 0..8 {
        if cube_values[i] < isolevel {
            cube_index |= 1 << i;
        }
    }

    // No triangles if all corners are on the same side
    if EDGE_TABLE[cube_index] == 0 {
        return;
    }

    // Find intersection points on edges
    let mut edge_vertices = [0usize; 12];

    for edge in 0..12 {
        if (EDGE_TABLE[cube_index] & (1 << edge)) != 0 {
            let (corner1, corner2) = EDGE_CORNERS[edge];
            let pos1 = get_corner_position(cube_pos, corner1);
            let pos2 = get_corner_position(cube_pos, corner2);
            let val1 = cube_values[corner1];
            let val2 = cube_values[corner2];

            // Create cache key for this edge
            let cache_key = create_edge_cache_key(pos1, pos2, edge);

            // Check if we've already calculated this edge intersection
            if let Some(&vertex_index) = vertex_cache.get(&cache_key) {
                edge_vertices[edge] = vertex_index;
            } else {
                // Calculate intersection point using linear interpolation
                let intersection = interpolate_edge(pos1, pos2, val1, val2, isolevel);
                let vertex_index = mesh.vertices.len();
                mesh.vertices.push(intersection);
                vertex_cache.insert(cache_key, vertex_index);
                edge_vertices[edge] = vertex_index;
            }
        }
    }

    // Generate triangles based on the lookup table
    let mut i = 0;
    while i < 16 && TRI_TABLE[cube_index][i] != -1 {
        let triangle = [
            edge_vertices[TRI_TABLE[cube_index][i] as usize],
            edge_vertices[TRI_TABLE[cube_index][i + 1] as usize],
            edge_vertices[TRI_TABLE[cube_index][i + 2] as usize],
        ];
        mesh.triangles.push(triangle);
        i += 3;
    }
}

fn get_point_value(
    points: &[f32],
    dimensions: (usize, usize, usize),
    pos: (usize, usize, usize),
) -> f32 {
    let (width, height, _depth) = dimensions;
    let (x, y, z) = pos;
    let index = z * height * width + y * width + x;
    points[index]
}

fn get_corner_position(cube_pos: (usize, usize, usize), corner: usize) -> [f32; 3] {
    let (x, y, z) = cube_pos;
    let (dx, dy, dz) = CORNER_OFFSETS[corner];
    [(x as f32) + dx, (y as f32) + dy, (z as f32) + dz]
}

fn interpolate_edge(
    pos1: [f32; 3],
    pos2: [f32; 3],
    val1: f32,
    val2: f32,
    isolevel: f32,
) -> [f32; 3] {
    if (val1 - val2).abs() < 1e-6 {
        return pos1; // Avoid division by zero
    }

    let t = (isolevel - val1) / (val2 - val1);
    [
        pos1[0] + t * (pos2[0] - pos1[0]),
        pos1[1] + t * (pos2[1] - pos1[1]),
        pos1[2] + t * (pos2[2] - pos1[2]),
    ]
}

fn create_edge_cache_key(
    pos1: [f32; 3],
    pos2: [f32; 3],
    edge: usize,
) -> (usize, usize, usize, usize) {
    // Create a consistent key regardless of edge direction
    let (p1, p2) = if pos1[0] < pos2[0]
        || (pos1[0] == pos2[0] && pos1[1] < pos2[1])
        || (pos1[0] == pos2[0] && pos1[1] == pos2[1] && pos1[2] < pos2[2])
    {
        (pos1, pos2)
    } else {
        (pos2, pos1)
    };

    // Convert to integer coordinates for hashing
    (
        (p1[0] as usize) * 1000 + (p1[1] as usize) * 100 + (p1[2] as usize) * 10,
        (p2[0] as usize) * 1000 + (p2[1] as usize) * 100 + (p2[2] as usize) * 10,
        edge,
        0, // Padding to make it a 4-tuple
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_cube() {
        // Simple 2x2x2 grid with one corner inside
        let points = vec![
            0.0, 0.5, // z=0, y=0 row
            0.5, 0.5, // z=0, y=1 row
            0.5, 0.5, // z=1, y=0 row
            0.5, 0.5, // z=1, y=1 row
        ];

        let mesh = marching_cubes(&points, (2, 2, 2), 0.25);

        // Should generate some triangles since we have values both above and below isolevel
        assert!(!mesh.triangles.is_empty());
        assert!(!mesh.vertices.is_empty());
    }
}
