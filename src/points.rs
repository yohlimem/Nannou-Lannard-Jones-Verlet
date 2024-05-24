
use nannou::{lyon::geom::Angle, prelude::*};


#[derive(Debug, Clone, Copy)]
pub struct Points {
    pub simulate: bool,
    pub pos: Vec2,
    pub last_pos: Vec2,
    pub v: Vec2,
    pub a: Vec2,
    pub mass: f32,
    pub charge: f32,
    pub potential: f32,
    pub plate_radius: f32
}

impl Points {
    //pub const SIGMA:f32 = 10.0;
    //pub const EPSILON:f32 = 10.0;
    const DT:f32 = 0.01;
    pub fn new(x0: Vec2, v0: Vec2, a: Vec2, plate_radius: f32, charge: f32, simulate: bool) -> Self {
        // println!("{:?}", Self::SIGMA);
        
        if !simulate {
            return Points {
                simulate,
                pos: x0,
                v: v0,
                a,
                mass: -1.0,
                plate_radius,
                charge,
                potential: 0.0,
                last_pos: x0, 
            }
        }
        Points {
            simulate,
            pos: x0 + v0*Self::DT+a*Self::DT*Self::DT,
            v: v0,
            a,
            mass: 1.0,
            plate_radius,
            charge,
            potential: 0.0,
            last_pos: x0, 
        }
    }

    pub fn random_init(plate_radius: f32, charge: f32, temperature: f32) -> Self {
        let angle = random_f32()*2.0*PI;
        let x = vec2(
            angle.cos(),
            angle.sin(),
        ) * random_f32()*plate_radius;
        // T = 0.5mv^2
        // v^2 = T/0.5m
        let angle = random_f32()*2.0*PI;
        let v = (temperature/0.5).sqrt()
        * vec2(
            angle.cos(),
            angle.sin(),
        );
        let a = vec2(0.0, 0.0);
        // println!("vel: {}", v);
        Points {
            simulate: true,
            pos: x + v*Self::DT+a*Self::DT*Self::DT,
            v,
            a,
            mass: 1.0,
            plate_radius,
            charge,
            last_pos: x,
            potential: 0.0,
        }
    }

    pub fn solver(&mut self, dt: f32) {
        // println!("===========================");
        // println!("pos: {}", self.pos);
        // println!("last_pos: {}", self.last_pos);
        if self.a.length() > 1000.0 {
            self.a = self.a.normalize() * 100.0;
        }
        self.v = (self.pos - self.last_pos) / 1.0001;
        self.last_pos = self.pos;
        self.pos += self.v + self.a * dt * dt;
        self.a = vec2(0.0, 0.0);
    }


    pub fn force(r: f32, epsilon: f32, sigma: f32) -> f32 {

        let r2 = r;
        let rp2 = r*r;
        let rp4 = rp2*rp2;
        let rp8 = rp4*rp4;
        let rp13 = rp8*rp4*r;
        let rp7 = rp4*rp2*r;

        let sigma2 = sigma*sigma;
        let sigma4 = sigma2*sigma2;
        let sigma6 = sigma4*sigma2;
        let sigma12 = sigma6*sigma6;

        let f = -2.0 * epsilon * ((12.0 * (sigma12) / rp13) - (6.0 * (sigma6) / rp7));
        f
    }

    fn charge_force(&self, p: &Points, r: f32) -> f32 {
        let k = 10.0 * 1000.0;
        (k * -self.charge * p.charge) / (r*r)
    }

    pub fn dist_to(&self, p: &Points) -> f32 {
        let dist = self.pos.distance(p.pos);
        if dist == 0.0 {
            return 0.01;
        }
        dist
    }
    pub fn dist_to_squared(&self, p: &Points) -> f32 {
        self.pos.distance_squared(p.pos)
    }

    pub fn r_vector(&self, p: &Points) -> Vec2 {
        let r = (p.pos - self.pos).normalize();
        if r.is_nan() {
            return vec2(0.0, 0.0);
        }
        r
        
    }

    pub fn step(&mut self, p_l: &[Points], epsilon: f32, sigma: f32) -> f32 {
        if !self.simulate {return 0.0;}
        for p in p_l {
            if p.pos == self.pos {
                continue;
            }
            let r = self.dist_to(p);
            let force = Self::force(r, epsilon, sigma);
            let r_vec = self.r_vector(p);
            let new_a = (force / self.mass) * r_vec + (self.charge_force(p, r) / self.mass) * r_vec;
            self.a += new_a;
            if self.pos.length() > self.plate_radius{
                self.a -= self.pos * 5.0;
            }
        }
        let total_force = self.a.length() * self.mass;
        self.solver(Self::DT);
        return total_force;
    }
}
