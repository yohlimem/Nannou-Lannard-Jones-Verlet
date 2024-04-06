
use nannou::prelude::*;


#[derive(Debug, Clone, Copy)]
pub struct Points {
    pub pos: Vec2,
    pub v: Vec2,
    pub a: Vec2,
    pub mass: f32,
    pub charge: f32,
    pub plate_size: (f32, f32)
}

impl Points {
    const SIGMA:f32 = 10.0;
    const EPSILON:f32 = 100.0;
    pub fn new(x0: Vec2, v0: Vec2, a: Vec2, plate_size: (f32, f32), charge: f32) -> Self {
        Points {
            pos: x0,
            v: v0,
            a,
            mass: 1.0,
            plate_size,
            charge,
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
        }
    }

    pub fn solver(&mut self, dt: f32) {
        self.v += self.a * dt;
        self.pos += self.v * dt;
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

    fn force(&self, r: f32) -> f32 {

        let mut r2 = 0.0;
        if r > Self::SIGMA*1.07{
            r2 = r
        }
        else if r < Self::SIGMA*0.1 {
            r2 = Self::SIGMA*0.9
        }
        else if r < Self::SIGMA*1.07 {
            r2 = Self::SIGMA*1.07
        }

        -4.0 * Self::EPSILON * ((12.0 * (Self::SIGMA.powi(12)) / r2.powi(13)) - (6.0 * (Self::SIGMA.powi(6)) / r2.powi(7)))
    }
    fn charge_force(&self, p: &Points, r: f32) -> f32 {
        let k = 30.0*Self::EPSILON;
        if r < 1.0 {
            return  (k * -self.charge * p.charge) / 1.0;
        }
        (k * -self.charge * p.charge) / (r*r)
    }

    pub fn dist_to(&self, p: &Points) -> f32 {
        self.pos.distance(p.pos)
    }

    pub fn r_vector(&self, p: &Points) -> Vec2 {
        (p.pos - self.pos).normalize()
    }

    pub fn step(&mut self, p_l: &[Points]) {
        for p in p_l {
            if p.pos == self.pos {
                continue;
            }
            let r = self.dist_to(p);
            if r > 300.0 {
                continue;
            }
            let r_vec = self.r_vector(p);
            let charge_force = self.charge_force(p, r);
            let force = self.force(r);
            let new_a = (self.force(r) / self.mass) * r_vec + (self.charge_force(p, r) / self.mass) * r_vec;
            self.a += new_a;
            // if force.abs() > 1000.0 {
            //     println!("force: {}", force);
            // }
            // if charge_force.abs() > 1000.0{
            //     println!("charge: {}", self.charge)
            // }
            // if force + charge_force > 1000.0 {
            //     println!("force + charge: {}", force + charge_force);
            // }
            if new_a.length() > 1000.0 {
                println!("r: {}", r);
                println!("force: {}", force);
                println!("charge: {}", self.charge);

                println!("new a: {}", new_a);
            }
        }
        if self.a.length() > 1000.0 {
            println!("a: {}", self.a);
        }
        self.v *= 0.9992;
        // self.a -= self.a.normalize();

        self.solver(0.01);
    }
}
