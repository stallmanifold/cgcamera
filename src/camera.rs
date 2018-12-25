use cgmath as math;
use cgmath::{Vector3, Vector4, Matrix4, Quaternion};

use std::fmt;


#[derive(Copy, Clone, Debug)]
pub struct FrustumFov {
    pub near: f32,
    pub far: f32,
    pub fov: f32,
    pub aspect: f32,
}

impl FrustumFov {
    pub fn new(near: f32, far: f32, fov: f32, aspect: f32) -> FrustumFov {
        FrustumFov {
            near: near,
            far: far,
            fov: fov,
            aspect: aspect,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct CameraAttitude {
    pub origin: Vector3,
    pub forward: Vector4,
    pub right: Vector4,
    pub up: Vector4,
    pub axis: Quaternion,
}

impl CameraAttitude {
    pub fn new(
        origin: Vector3, forward: Vector4,
        right: Vector4, up: Vector4, rotation_axis: Vector3) -> CameraAttitude {

        let axis = Quaternion::from_sv(0.0, rotation_axis);
        CameraAttitude {
            origin: origin,
            forward: forward,
            right: right,
            up: up,
            axis: axis,
        }
    }
}

#[derive(Clone, Debug)]
pub struct PerspectiveFovCamera {
    // Camera parameters.
    pub near: f32,
    pub far: f32,
    pub fov: f32,
    pub aspect: f32,

    // Camera kinematics.
    pub origin: Vector3,
    pub forward: Vector4,
    pub right: Vector4,
    pub up: Vector4,
    pub axis: Quaternion,

    // Camera matrices.
    pub proj_mat: Matrix4,
    pub trans_mat: Matrix4,
    pub rot_mat: Matrix4,
    pub view_mat: Matrix4,
}

impl PerspectiveFovCamera {
    pub fn new(frustum: FrustumFov, attitude: CameraAttitude) -> PerspectiveFovCamera {
        let proj_mat = math::perspective((
            frustum.fov, frustum.aspect, frustum.near, frustum.far
        ));
        let trans_mat = Matrix4::from_translation(attitude.origin);
        let rot_mat = Matrix4::from(attitude.axis);
        let view_mat = rot_mat * trans_mat;

        PerspectiveFovCamera {
            near: frustum.near,
            far: frustum.far,
            fov: frustum.fov,
            aspect: frustum.aspect,

            origin: attitude.origin,
            forward: attitude.forward,
            right: attitude.right,
            up: attitude.up,
            axis: attitude.axis,

            proj_mat: proj_mat,
            trans_mat: trans_mat,
            rot_mat: rot_mat,
            view_mat: view_mat,
        }
    }
}

impl fmt::Display for PerspectiveFovCamera {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Camera Model:").unwrap();
        writeln!(f, "near: {}", self.near).unwrap();
        writeln!(f, "far: {}", self.far).unwrap();
        writeln!(f, "aspect: {}", self.aspect).unwrap();
        writeln!(f, "origin: {}", self.origin).unwrap();
        writeln!(f, "fwd: {}", self.forward).unwrap();
        writeln!(f, "rgt: {}", self.right).unwrap();
        writeln!(f, "up: {}", self.up).unwrap();
        writeln!(f, "axis: {}", self.axis).unwrap();
        writeln!(f, "proj_mat: {}", self.proj_mat).unwrap();
        writeln!(f, "trans_mat: {}", self.trans_mat).unwrap();
        writeln!(f, "rot_mat: {}", self.rot_mat).unwrap();
        writeln!(f, "view_mat: {}", self.view_mat)
    }
}


#[derive(Copy, Clone, Debug)]
pub struct Frustum {
    pub left: f32,
    pub right: f32,
    pub bottom: f32,
    pub top: f32,
    pub near: f32,
    pub far: f32,
}

impl Frustum {
    pub fn new(
        left: f32, right: f32,
        bottom: f32, top: f32, near: f32, far: f32) -> Frustum {

        Frustum {
            left: left,
            right: right,
            bottom: bottom,
            top: top,
            near: near,
            far: far,
        }
    }
}


#[derive(Clone, Debug)]
pub struct PerspectiveCamera {
    // Camera frustum parameters.
    pub left: f32,
    pub right: f32,
    pub bottom: f32,
    pub top: f32,
    pub near: f32,
    pub far: f32,

    // Camera kinematics.
    pub origin: Vector3,
    pub forward: Vector4,
    pub rgt: Vector4,
    pub up: Vector4,
    pub axis: Quaternion,

    // Camera matrices.
    pub proj_mat: Matrix4,
    pub trans_mat: Matrix4,
    pub rot_mat: Matrix4,
    pub view_mat: Matrix4,
}

impl PerspectiveCamera {
    pub fn new(frustum: Frustum, attitude: CameraAttitude) -> PerspectiveCamera {
        let proj_mat = math::frustum((
            frustum.left, frustum.right,
            frustum.bottom, frustum.top, frustum.near, frustum.far
        ));
        let trans_mat = Matrix4::from_translation(attitude.origin);
        let rot_mat = Matrix4::from(attitude.axis);
        let view_mat = rot_mat * trans_mat;

        PerspectiveCamera {
            left: frustum.left,
            right: frustum.right,
            bottom: frustum.bottom,
            top: frustum.top,
            near: frustum.near,
            far: frustum.far,

            origin: attitude.origin,
            forward: attitude.forward,
            rgt: attitude.right,
            up: attitude.up,
            axis: attitude.axis,

            proj_mat: proj_mat,
            trans_mat: trans_mat,
            rot_mat: rot_mat,
            view_mat: view_mat,
        }
    }
}

impl fmt::Display for PerspectiveCamera {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Camera Model:").unwrap();
        writeln!(f, "left: {}", self.left).unwrap();
        writeln!(f, "right: {}", self.right).unwrap();
        writeln!(f, "bottom: {}", self.bottom).unwrap();
        writeln!(f, "top: {}", self.top).unwrap();
        writeln!(f, "near: {}", self.near).unwrap();
        writeln!(f, "far: {}", self.far).unwrap();
        writeln!(f, "origin: {}", self.origin).unwrap();
        writeln!(f, "fwd: {}", self.forward).unwrap();
        writeln!(f, "rgt: {}", self.right).unwrap();
        writeln!(f, "up: {}", self.up).unwrap();
        writeln!(f, "axis: {}", self.axis).unwrap();
        writeln!(f, "proj_mat: {}", self.proj_mat).unwrap();
        writeln!(f, "trans_mat: {}", self.trans_mat).unwrap();
        writeln!(f, "rot_mat: {}", self.rot_mat).unwrap();
        writeln!(f, "view_mat: {}", self.view_mat)
    }
}

#[derive(Clone, Debug)]
pub struct OrthographicCamera {
    // Camera frustum parameters.
    pub left: f32,
    pub right: f32,
    pub bottom: f32,
    pub top: f32,
    pub near: f32,
    pub far: f32,

    // Camera kinematics.
    pub origin: Vector3,
    pub forward: Vector4,
    pub rgt: Vector4,
    pub up: Vector4,
    pub axis: Quaternion,

    // Camera matrices.
    pub proj_mat: Matrix4,
    pub trans_mat: Matrix4,
    pub rot_mat: Matrix4,
    pub view_mat: Matrix4,
}

impl OrthographicCamera {
    pub fn new(frustum: Frustum, attitude: CameraAttitude) -> OrthographicCamera {
        let proj_mat = math::ortho((
            frustum.left, frustum.right,
            frustum.bottom, frustum.top, frustum.near, frustum.far
        ));
        let trans_mat = Matrix4::from_translation(attitude.origin);
        let rot_mat = Matrix4::from(attitude.axis);
        let view_mat = rot_mat * trans_mat;

        OrthographicCamera {
            left: frustum.left,
            right: frustum.right,
            bottom: frustum.bottom,
            top: frustum.top,
            near: frustum.near,
            far: frustum.far,

            origin: attitude.origin,
            forward: attitude.forward,
            rgt: attitude.right,
            up: attitude.up,
            axis: attitude.axis,

            proj_mat: proj_mat,
            trans_mat: trans_mat,
            rot_mat: rot_mat,
            view_mat: view_mat,
        }
    }
}

impl fmt::Display for OrthographicCamera {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Camera Model:").unwrap();
        writeln!(f, "left: {}", self.left).unwrap();
        writeln!(f, "right: {}", self.right).unwrap();
        writeln!(f, "bottom: {}", self.bottom).unwrap();
        writeln!(f, "top: {}", self.top).unwrap();
        writeln!(f, "near: {}", self.near).unwrap();
        writeln!(f, "far: {}", self.far).unwrap();
        writeln!(f, "origin: {}", self.origin).unwrap();
        writeln!(f, "fwd: {}", self.forward).unwrap();
        writeln!(f, "rgt: {}", self.rgt).unwrap();
        writeln!(f, "up: {}", self.up).unwrap();
        writeln!(f, "axis: {}", self.axis).unwrap();
        writeln!(f, "proj_mat: {}", self.proj_mat).unwrap();
        writeln!(f, "trans_mat: {}", self.trans_mat).unwrap();
        writeln!(f, "rot_mat: {}", self.rot_mat).unwrap();
        writeln!(f, "view_mat: {}", self.view_mat)
    }
}


#[cfg(test)]
mod orthographic_camera_tests {
    use super::{CameraAttitude, Frustum, OrthographicCamera};
    use cgmath::{Vector3, Vector4, Matrix4};

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
    fn test_unit_viewing_volume_projection_matrix_should_be_mirroring_matrix() {
        let left = -1.0;
        let right = 1.0;
        let bottom = -1.0;
        let top = 1.0;
        let near = 1.0;
        let far = -1.0;
        let frustum = Frustum::new(left, right, bottom, top, near, far);

        let origin = cgmath::vec3((0.0, 0.0, 5.0));
        let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
        let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
        let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
        let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
        let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
        let camera = OrthographicCamera::new(frustum, attitude);

        let mirror_mat = Matrix4::new(
            1.0, 0.0,  0.0, 0.0,
            0.0, 1.0,  0.0, 0.0,
            0.0, 0.0, -1.0, 0.0,
            0.0, 0.0,  0.0, 1.0
        );

        assert_eq!(camera.proj_mat, mirror_mat);
    }

    #[test]
    fn test_orthographic_camera_should_map_points_to_canonical_view_volume() {
        let left = -4.0;
        let right = 4.0;
        let bottom = -4.0;
        let top = 4.0;
        let near = 4.0;
        let far = -4.0;
        let frustum = Frustum::new(left, right, bottom, top, near, far);

        let origin = cgmath::vec3((0.0, 0.0, 5.0));
        let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
        let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
        let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
        let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
        let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
        let camera = OrthographicCamera::new(frustum, attitude);

        let p_wor = cgmath::vec4((-3.0, 3.0, 3.0, 1.0));
        let p_cam = camera.view_mat * p_wor;
        let p_cvv = camera.proj_mat * p_cam;
        let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
        let cvv = AxisAlignedBoundingBox::unit_aabb();

        assert!(cvv.contains(p_cvv), "p_wor = {}; p_cvv = {}", p_wor, p_cvv);
    }
}

