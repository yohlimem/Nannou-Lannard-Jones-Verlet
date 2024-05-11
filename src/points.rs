
use nannou::prelude::*;


#[derive(Debug, Clone, Copy)]
pub struct Points {
    pub pos: Vec2,
    pub last_pos: Vec2,
    pub v: Vec2,
    pub a: Vec2,
    pub mass: f32,
    pub charge: f32,
    pub potential: f32,
    pub plate_size: (f32, f32)
}

impl Points {
    //pub const SIGMA:f32 = 10.0;
    //pub const EPSILON:f32 = 10.0;
    pub fn new(x0: Vec2, v0: Vec2, a: Vec2, plate_size: (f32, f32), charge: f32) -> Self {
        // println!("{:?}", Self::SIGMA);
        Points {
            pos: x0,
            v: v0,
            a,
            mass: 1.0,
            plate_size,
            charge,
            potential: 0.0,
            last_pos: x0,
        }
    }

    pub fn random_init(plate_size: Vec2, charge: f32) -> Self {
        let x = vec2(
            random_range(0.0, plate_size.x),
            random_range(0.0, plate_size.y),
        );
        let v = vec2(0.0, 0.0);
        let a = vec2(0.0, 0.0);

        Points {
            pos: x,
            v,
            a,
            mass: 1.0,
            plate_size: (plate_size.x, plate_size.y),
            charge,
            last_pos: x,
            potential: 0.0,
        }
    }

    pub fn solver(&mut self, dt: f32) {
        // self.v += self.a * dt;
        self.pos += self.pos - self.last_pos + self.a * dt * dt;
        self.last_pos = self.pos;
        // println!("{:?}", self.pos);
        if self.pos.x > self.plate_size.0 || self.pos.x < 0.0 {
            self.v.x *= -1.0;
            // self.pos += self.v * dt;
        }

        if self.pos.y > self.plate_size.1 || self.pos.y < 0.0 {
            self.v.y *= -1.0;
            // self.pos += self.v * dt;
        }
        self.a = vec2(0.0, 0.0);
    }

    pub fn force(r: f32, epsilon: f32, sigma: f32) -> f32 {

        let r2 = r;
        // if r > Self::SIGMA*1.07{
        //     r2 = r
        // }
        // else if r < 0.1{
        //     r2 = Self::SIGMA*100.0;
        // }
        // else if r < Self::SIGMA*0.1 {
        //     r2 = Self::SIGMA*2.0
        // }
        // else if r < Self::SIGMA*1.07 {
        //     r2 = Self::SIGMA*1.07
        // }

        if r2.powi(13) <= 0.001 {
            println!("die");
        }


        let f = -2.0 * epsilon * ((12.0 * (sigma.powi(12)) / r2.powi(13)) - (6.0 * (sigma.powi(6)) / r2.powi(7)));
        if f.abs() > 1000.0 {
            println!("f: {}", f);
            println!("A1: {}", (12.0 * (sigma.powi(12))) / r2.powi(13));
            println!("A2: {}", (6.0 * (sigma.powi(6)) / r2.powi(7)));
        }
        f
    }

    fn u(&self, r: f32) -> f32 {
        // let mut r2 = r;
        // if r > Self::SIGMA*1.07{
        //     r2 = r
        // }
        // else if r < Self::SIGMA*0.1 {
        //     r2 = Self::SIGMA*2.0
        // }
        // else if r < Self::SIGMA*1.07 {
        //     r2 = Self::SIGMA*1.07
        // }
        //2.0 * Self::EPSILON * ((Self::SIGMA.powi(12) / r.powi(12)) - (Self::SIGMA.powi(6) / r.powi(6)))
        0.0
    }

    fn charge_force(&self, p: &Points, r: f32) -> f32 {
        let k = 10.0 * 1000.0;
        // if r < 1.0 {
        //     return  (k * -self.charge * p.charge) / 1.0;
        // }
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
        for p in p_l {
            if p.pos == self.pos {
                continue;
            }
            let r = self.dist_to(p);
            if r > 100.0 {
                continue;
            }
            let force = Self::force(r, epsilon, sigma);
            let r_vec = self.r_vector(p);
            let new_a = (force / self.mass) * r_vec + (self.charge_force(p, r) / self.mass) * r_vec;
            // self.potential += self.u(r);
            // print!("{}", self.potential);
            if r_vec.is_nan() {
                println!("r_vec: {}", r_vec);
            }
            if force.is_nan() {
                println!("force: nan");
            }
            // if self.charge_force(p, r).is_nan() {
            //     println!("charge: {}", self.charge)
            // }
            // if new_a.length() > 1000.0 {
                
            // }
            self.a += new_a;
        }
        let total_force = self.a.length() * self.mass;
        self.solver(0.1);
        return total_force;
    }
}
