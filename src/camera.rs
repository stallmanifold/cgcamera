use cgmath as math;
use cgmath::{Vector3, Vector4, Matrix4, Quaternion};

use std::fmt;


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

#[derive(Clone, Debug)]
pub struct PerspectiveFovCamera {
    // Camera parameters.
    pub near: f32,
    pub far: f32,
    pub fov: f32,
    pub aspect: f32,

    // Camera kinematics.
    pub origin: Vector3,
    pub fwd: Vector4,
    pub rgt: Vector4,
    pub up: Vector4,
    pub axis: Quaternion,

    // Camera matrices.
    pub proj_mat: Matrix4,
    pub trans_mat: Matrix4,
    pub rot_mat: Matrix4,
    pub view_mat: Matrix4,
}

impl PerspectiveFovCamera {
    pub fn new(
        frustum: FrustumFov,
        cam_pos: Vector3,
        fwd: Vector4, rgt: Vector4, up: Vector4, axis: Vector3) -> PerspectiveFovCamera {

        let proj_mat = math::perspective(
            (frustum.fov, frustum.aspect, frustum.near, frustum.far)
        );
        let trans_mat = Matrix4::from_translation(cam_pos);
        let axis_quat = Quaternion::from_sv(0.0, axis);
        let rot_mat = Matrix4::from(axis_quat);
        let view_mat = rot_mat * trans_mat;

        PerspectiveFovCamera {
            near: frustum.near,
            far: frustum.far,
            fov: frustum.fov,
            aspect: frustum.aspect,

            origin: cam_pos,
            fwd: fwd,
            rgt: rgt,
            up: up,
            axis: axis_quat,

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
        writeln!(f, "fwd: {}", self.fwd).unwrap();
        writeln!(f, "rgt: {}", self.rgt).unwrap();
        writeln!(f, "up: {}", self.up).unwrap();
        writeln!(f, "axis: {}", self.axis).unwrap();
        writeln!(f, "proj_mat: {}", self.proj_mat).unwrap();
        writeln!(f, "trans_mat: {}", self.trans_mat).unwrap();
        writeln!(f, "rot_mat: {}", self.rot_mat).unwrap();
        writeln!(f, "view_mat: {}", self.view_mat)
    }
}


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
    pub fwd: Vector4,
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
    pub fn new(
        near: f32, far: f32, bottom: f32, top: f32, left: f32, right: f32,
        cam_pos: Vector3,
        fwd: Vector4, rgt: Vector4, up: Vector4, axis: Vector3) -> PerspectiveCamera {

        let proj_mat = math::frustum((left, right, bottom, top, near, far));
        let trans_mat = Matrix4::from_translation(cam_pos);
        let axis_quat = Quaternion::from_sv(0.0, axis);
        let rot_mat = Matrix4::from(axis_quat);
        let view_mat = rot_mat * trans_mat;

        PerspectiveCamera {
            left: left,
            right: right,
            bottom: bottom,
            top: top,
            near: near,
            far: far,

            origin: cam_pos,
            fwd: fwd,
            rgt: rgt,
            up: up,
            axis: axis_quat,

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
        writeln!(f, "fwd: {}", self.fwd).unwrap();
        writeln!(f, "rgt: {}", self.rgt).unwrap();
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
    pub fwd: Vector4,
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
    pub fn new(
        near: f32, far: f32, bottom: f32, top: f32, left: f32, right: f32,
        cam_pos: Vector3,
        fwd: Vector4, rgt: Vector4, up: Vector4, axis: Vector3) -> OrthographicCamera {

        let proj_mat = math::frustum((left, right, bottom, top, near, far));
        let trans_mat = Matrix4::from_translation(cam_pos);
        let axis_quat = Quaternion::from_sv(0.0, axis);
        let rot_mat = Matrix4::from(axis_quat);
        let view_mat = rot_mat * trans_mat;

        OrthographicCamera {
            left: left,
            right: right,
            bottom: bottom,
            top: top,
            near: near,
            far: far,

            origin: cam_pos,
            fwd: fwd,
            rgt: rgt,
            up: up,
            axis: axis_quat,

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
        writeln!(f, "fwd: {}", self.fwd).unwrap();
        writeln!(f, "rgt: {}", self.rgt).unwrap();
        writeln!(f, "up: {}", self.up).unwrap();
        writeln!(f, "axis: {}", self.axis).unwrap();
        writeln!(f, "proj_mat: {}", self.proj_mat).unwrap();
        writeln!(f, "trans_mat: {}", self.trans_mat).unwrap();
        writeln!(f, "rot_mat: {}", self.rot_mat).unwrap();
        writeln!(f, "view_mat: {}", self.view_mat)
    }
}
