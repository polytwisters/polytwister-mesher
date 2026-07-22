use std::default;
use serde::Deserialize;


#[derive(Deserialize, Clone, Copy, Debug)]
#[serde(deny_unknown_fields)]
pub struct CylinderMeshConfig {
    #[serde(default = "CylinderMeshConfig::default_half_length")]
    pub half_length: f64,
    #[serde(default = "CylinderMeshConfig::default_linear_segments")]
    pub linear_segments: usize,
    #[serde(default = "CylinderMeshConfig::default_radial_segments")]
    pub radial_segments: usize,

    #[serde(default = "CylinderMeshConfig::default_resolution")]
    pub resolution: f64,
}

impl CylinderMeshConfig {
    fn default_half_length() -> f64 { 2.0 }
    fn default_linear_segments() -> usize { 50usize }
    fn default_radial_segments() -> usize { 50usize }
    fn default_resolution() -> f64 { 0.05 }
}

impl Default for CylinderMeshConfig {
    fn default() -> Self {
        Self {
            half_length: Self::default_half_length(),
            linear_segments: Self::default_linear_segments(),
            radial_segments: Self::default_radial_segments(),
            resolution: Self::default_resolution(),
        }
    }
}

#[derive(Deserialize, Clone, Copy, Debug)]
#[serde(deny_unknown_fields)]
pub struct TorusMeshConfig {
    #[serde(default = "TorusMeshConfig::default_radius")]
    pub radius: f64,
    #[serde(default = "TorusMeshConfig::default_linear_segments")]
    pub linear_segments: usize,
    #[serde(default = "TorusMeshConfig::default_radial_segments")]
    pub radial_segments: usize,
}

impl TorusMeshConfig {
    fn default_radius() -> f64 { 0.01 }
    fn default_linear_segments() -> usize { 50usize }
    fn default_radial_segments() -> usize { 16usize }
}

impl Default for TorusMeshConfig {
    fn default() -> Self {
        Self {
            radius: Self::default_radius(),
            linear_segments: Self::default_linear_segments(),
            radial_segments: Self::default_radial_segments(),
        }
    }
}

#[derive(Deserialize, Clone, Copy, Debug)]
#[serde(deny_unknown_fields)]
pub struct RingMeshConfig {
    #[serde(default = "RingMeshConfig::default_radius")]
    pub radius: f64,
    #[serde(default = "RingMeshConfig::default_longitudes")]
    pub longitudes: usize,
    #[serde(default = "RingMeshConfig::default_latitudes")]
    pub latitudes: usize,
}

impl RingMeshConfig {
    fn default_radius() -> f64 { 0.02 }
    fn default_longitudes() -> usize { 16 }
    fn default_latitudes() -> usize { 16 }
}

impl Default for RingMeshConfig {
    fn default() -> Self {
        Self {
            radius: Self::default_radius(),
            longitudes: Self::default_longitudes(),
            latitudes: Self::default_latitudes(),
        }
    }
}

#[derive(Deserialize, Clone, Copy, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct Config {
    #[serde(default)]
    pub rings: RingMeshConfig,
    #[serde(default)]
    pub strips: TorusMeshConfig,
    #[serde(default)]
    pub twisters: CylinderMeshConfig,
}
