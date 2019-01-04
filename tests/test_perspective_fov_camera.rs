use cgmath::Vector3;
use cgcamera::{CameraAttitude, FrustumFov, PerspectiveFovCamera};


struct AxisAlignedBoundingBox {
    left: f32,
    right: f32,
    bottom: f32,
    top: f32,
    near: f32,
    far: f32,
}

impl AxisAlignedBoundingBox {
    fn new(left: f32, right: f32, bottom: f32, top: f32, near: f32, far: f32) -> Self {
        Self { left, right, bottom, top, near, far }
    }

    fn contains(&self, point: Vector3) -> bool {
        (point.x >= self.left) && (point.x <= self.right) &&
            (point.y >= self.bottom) && (point.y <= self.top) &&
            (point.z >= self.far && point.z <= self.near)
    }

    fn unit_aabb() -> Self {
        Self::new(-1.0, 1.0, -1.0, 1.0, 1.0, -1.0)
    }
}

#[test]
fn test_perspective_camera_should_map_points_inside_frustum_inside_canonical_view_volume() {
    let near = 0.1;
    let far = 100.0;
    let fovy = 67.0; // degrees
    let aspect = 720.0 / 480.0;
    let frustum_fov = FrustumFov::new(near, far, fovy, aspect);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
    let camera = PerspectiveFovCamera::new(frustum_fov, attitude);

    let p_cam = cgmath::vec4((1.0, 1.0, 1.0, 1.0));
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(cvv.contains(p_cvv), "p_cam = {}; p_cvv = {}", p_cam, p_cvv);
}

#[test]
fn test_perspective_camera_should_map_points_outside_frustum_outside_canonical_view_volume() {
    let near = 0.1;
    let far = 100.0;
    let fovy = 67.0; // degrees
    let fovy_rad = fovy * cgmath::ONE_DEG_IN_RAD;
    let aspect = 720.0 / 480.0;
    let frustum_fov = FrustumFov::new(near, far, fovy, aspect);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
    let camera = PerspectiveFovCamera::new(frustum_fov, attitude);

    let p_cam = cgmath::vec4((300.0, 300.0, 300.0, 1.0));
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(!cvv.contains(p_cvv), "p_cam = {}; p_cvv = {}", p_cam, p_cvv);
}
