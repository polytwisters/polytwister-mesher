use serde::Deserialize;


#[derive(Deserialize, Clone, Copy, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct CylinderMeshConfig {
    #[serde(default = "CylinderMeshConfig::default_half_length")]
    pub half_length: f64,
    #[serde(default = "CylinderMeshConfig::default_linear_segments")]
    pub linear_segments: usize,
    #[serde(default = "CylinderMeshConfig::default_radial_segments")]
    pub radial_segments: usize,
}

impl CylinderMeshConfig {
    fn default_half_length() -> f64 { 2.0 }
    fn default_linear_segments() -> usize { 50usize }
    fn default_radial_segments() -> usize { 50usize }
}

#[derive(Deserialize, Clone, Copy, Debug, Default)]
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

#[derive(Deserialize, Clone, Copy, Debug, Default)]
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
