use serde::Deserialize;


#[derive(Deserialize)]
#[serde(rename_all="camelCase", default)]
#[derive(Clone, Copy, Debug)]
pub struct CylinderMeshConfig {
    pub half_length: f64,
    pub linear_segments: usize,
    pub radial_segments: usize,
}

impl Default for CylinderMeshConfig {
    fn default() -> Self {
        CylinderMeshConfig {
            half_length: 5.0,
            linear_segments: 100,
            radial_segments: 100,
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all="camelCase", default)]
#[derive(Clone, Copy, Debug)]
pub struct TorusMeshConfig {
    pub thickness: f64,
    pub linear_segments: usize,
    pub radial_segments: usize,
}

impl Default for TorusMeshConfig {
    fn default() -> Self {
        TorusMeshConfig {
            thickness: 0.01,
            radial_segments: 16,
            linear_segments: 128,
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all="camelCase", default)]
#[derive(Clone, Copy, Debug)]
pub struct RingMeshConfig {
    pub radius: f64,
    pub segments: usize,
    pub rings: usize,
}

impl Default for RingMeshConfig {
    fn default() -> Self {
        RingMeshConfig {
            radius: 0.02,
            segments: 16,
            rings: 32, 
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all="camelCase", default)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Config {
    pub rings: RingMeshConfig,
    pub strips: TorusMeshConfig,
    pub twisters: CylinderMeshConfig,
}
