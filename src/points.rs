
use nannou::{lyon::geom::Angle, prelude::*};
use std::collections::vec_deque::VecDeque;


#[derive(Debug, Clone)]
pub struct Points {
    pub simulate: bool,
    pub pos: Vec2,
    pub last_pos: Vec2,
    pub v: Vec2,
    pub a: Vec2,
    pub mass: f32,
    pub charge: f32,
    pub potential: f32,
    pub plate_radius: f32,
    pub temperature: f32,
    pub got_good: bool,
    pub outside: bool,
    pub last_ten_velocities: VecDeque<Vec2>,

}

impl Points {
    //pub const SIGMA:f32 = 10.0;
    //pub const EPSILON:f32 = 10.0;
    pub const DT:f32 = 0.01;
    pub fn new(x0: Vec2, v0: Vec2, a: Vec2, plate_radius: f32, charge: f32, simulate: bool) -> Self { // create a new point
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
                temperature: v0.length()*v0.length()*0.5,
                got_good: false,
                last_ten_velocities: VecDeque::with_capacity(10),
                outside: false,
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
            temperature: v0.length()*v0.length()*0.5,
            got_good: false,
            last_ten_velocities: VecDeque::with_capacity(10),
            outside: false,

        }
    }

    pub fn random_init(plate_radius: f32, charge: f32, temperature: f32) -> Self { // create a random new point using a radius
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
            pos: x + v*Self::DT+a*Self::DT*Self::DT*0.5,
            v,
            temperature: v.length()*v.length()*0.5,

            a,
            mass: 1.0,
            plate_radius,
            charge,
            last_pos: x,
            potential: 0.0,
            got_good: false,
            last_ten_velocities: VecDeque::with_capacity(10),
            outside: false,
        }
    }

    pub fn solver(&mut self, dt: f32, temperature_depletion: f32) { // move the points according to an integration method
        // println!("===========================");
        // println!("pos: {}", self.pos);
        // println!("last_pos: {}", self.last_pos);
        // println!("temp: {}", self.temperature);
        // println!("vel: {}", self.v.length());
        if self.a.length() > 1000.0 {
            self.a = self.a.normalize() * 1000.0;
        }
        self.v = (self.pos - self.last_pos) * temperature_depletion;
        if self.outside {
            // println!("vel: {}", self.v);
            let normal = -self.pos.normalize();
            self.v -= 2.0 * self.v.dot(normal) / normal.dot(normal) * normal;
            self.pos = (self.plate_radius - 0.1) * self.pos.normalize();
        }
        self.last_pos = self.pos;
        self.pos += self.v + self.a * dt * dt;
        self.a = vec2(0.0, 0.0);
        self.temperature = self.mass*self.v.length()*self.v.length()*0.5;
        // self.v = Vec2::ZERO;
    }


    fn force(r: f32, epsilon: f32, sigma: f32) -> f32 { // force between two points using the Lannard Jones potential

        let rp2 = r*r;
        let rp4 = rp2*rp2;
        let rp8 = rp4*rp4;
        let rp13 = rp8*rp4*r;
        let rp7 = rp4*rp2*r;

        let sigma2 = sigma*sigma;
        let sigma4 = sigma2*sigma2;
        let sigma6 = sigma4*sigma2;
        let sigma12 = sigma6*sigma6;

        let f = -4.0 * epsilon * ((12.0 * (sigma12) / rp13) - (6.0 * (sigma6) / rp7));
        f
    }

    fn force_with_r_squared(r_squared: f32, epsilon: f32, sigma: f32) -> f32 { // use r squared to avoid calculating the square root (which didnt happen lol)
        let r = r_squared.sqrt();
        let rp4 = r_squared*r_squared;
        let rp8 = rp4*rp4;
        let rp13 = rp8*rp4*r;
        let rp7 = rp4*r_squared*r;

        let sigma2 = sigma*sigma;
        let sigma4 = sigma2*sigma2;
        let sigma6 = sigma4*sigma2;
        let sigma12 = sigma6*sigma6;

        let f = -4.0 * epsilon * ((12.0 * (sigma12) / rp13) - (6.0 * (sigma6) / rp7));
        f
    }


    fn charge_force(&self, p: &Points, r: f32) -> f32 { // calculate the force between two points using the coulomb potential which is the magnetic force between two charges
        let k = 10.0 * 1000.0;
        (k * -self.charge * p.charge) / (r*r)
    }

    fn charge_force_r_squared(&self, p: &Points, r_squared: f32) -> f32 { // calculate the force between two points using the coulomb potential which is the magnetic force between two charges using r squared to avoid calculating the square root
        let k = 10.0 * 1000.0;
        (k * -self.charge * p.charge) / r_squared
    }

    fn dist_to(&self, p: &Points) -> f32 { // distance between two points
        let dist = self.pos.distance(p.pos);
        if dist == 0.0 {
            return 0.01;
        }
        dist
    }
    fn dist_to_squared(&self, p: &Points) -> f32 { // distance squared between two points
        self.pos.distance_squared(p.pos)
    }

    fn r_vector(&self, p: &Points) -> Vec2 { // vector between two points
        let r = p.pos - self.pos;
        if r.is_nan() {
            return vec2(0.0, 0.0);
        }
        r
        
    }
    

    pub fn step(&mut self, p_l: &[Points], epsilon: f32, sigma: f32 , temperature_depletion: f32) { // do a step of the simulation
        if !self.simulate {return;}
        for p in p_l {
            if p.pos == self.pos {
                continue;
            }
            let r_vec = self.r_vector(p); // calculate the vector between two points
            let r = r_vec.length_squared(); // LENGTH SQUARED!!!

            let force = Self::force_with_r_squared(r, epsilon, sigma); // use the Lannard Jones potential to calculate the force between two points its
            let new_a = (force + self.charge_force_r_squared(p, r)) * r_vec.normalize() / self.mass; // calculate the acceleration of the point
            self.a += new_a; // add the acceleration to the point
            
            self.outside = self.pos.length() > self.plate_radius
        }


        self.solver(Self::DT, temperature_depletion); // use the solver to move the points
    }

    pub fn test_crystalized(&mut self){
        self.last_ten_velocities.push_front(self.v);
        if self.last_ten_velocities.len() > 20 {
            self.last_ten_velocities.pop_back();
        }
        let mut sum_temp = Vec2::ZERO;
        for v in self.last_ten_velocities.iter() {
            sum_temp += *v;
        }
        self.got_good = (sum_temp.length_squared()/20.0) < 0.05;
    }
}
