use nannou::prelude::*;

use crate::Model;

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
        let epsilon = 100.0; // replace with actual value
        let sigma = 10.0; // replace with actual value

        let mut r2 = 0.0;
        if r > sigma*1.07{
            r2 = r
        }
        else{
            r2 = sigma*1.07
        }

        -4.0 * epsilon * ((12.0 * (sigma.powi(12)) / r2.powi(13)) - (6.0 * (sigma.powi(6)) / r2.powi(7)))
    }
    fn charge_force(&self, p: &Points, r: f32) -> f32 {
        let k = 1000.0;
        if r < 2.0 {
            return  k * -self.charge * p.charge / 1.0;
        }
        k * -self.charge * p.charge / (r*r)
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
            if r > 100.0 {
                continue;
            }
            self.a += (self.force(r) / self.mass) * self.r_vector(p);
            self.a += (self.charge_force(p, r) / self.mass) * self.r_vector(p);
        }
        self.v *= 0.9992;
        // self.a -= self.a.normalize();

        self.solver(0.01);
    }
}
