use egui_plot::{Line, Plot, PlotPoints};
use nannou::{color::rgb::Rgb, draw::theme::Color, prelude::*, text::rt::Point};
use nannou_egui::{
    self,
    egui::{self, Response},
    Egui,
};
use rand::Rng;
mod points;
use points::Points;
/*sigma = 10 # nm
epsilon = 500 # J
K = 20 # temporary
plate_size = (500, 500) # nm^2
p_l = [(100,100), (60, 35), (20,200), (300, 300)] # list of all particles
         # example: [(0, 0), (30, 403)]
p_l_copy = copy.deepcopy(p_l)
 */
struct Model {
    // window: Window,
    egui: Egui,
    sigma: f32,             // nm
    epsilon: f32,           // J
    K: f32,                 // temporary
    plate_size: (f32, f32), // nm^2
    p_l: Vec<Points>,       // list of all particles
    p_l_copy: Vec<Points>,
    speed: f32,
    time: f32,
    done_simulation: bool,
    total_force: f32,
    seed_size: u32,
    atom_count: u32,
}

fn gen_atoms(model: &mut Model, n: u32, seed_size: u32) {
    let mut l: Vec<Points> = Vec::new();
    let mut rng = rand::thread_rng();

    let seed_x0 = model.plate_size.0 / 2.0
        + (0.0 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
    let seed_y0 = model.plate_size.1 / 2.0
        + (0.0 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
    let seed_x1 =
        model.plate_size.0 / 2.0 + (f32::sqrt(seed_size as f32) - 1.0) / 2.0 * 1.2 * model.sigma;
    let seed_y1 =
        model.plate_size.1 / 2.0 + (f32::sqrt(seed_size as f32) - 1.0) / 2.0 * 1.2 * model.sigma;
    for i in 0..n {
        let x = rng.gen_range(0.0..model.plate_size.0);
        let y = rng.gen_range(0.0..model.plate_size.1);
        if !(x < seed_x0 - 20.0 || x > seed_x1 + 20.0)
            && !(y < seed_y0 - 20.0 || y > seed_y1 + 20.0)
        {
            continue;
        }

        let p = Points::new(
            vec2(x, y),
            vec2(0.0, 0.0),
            vec2(0.0, 0.0),
            model.plate_size,
            i as f32 % 2.0 - 0.4,
        );

        l.push(p);
    }
    for x in 0..(f32::sqrt(seed_size as f32) as u32) {
        for y in 0..(f32::sqrt(seed_size as f32) as u32) {
            
            let new_x = model.plate_size.0 / 2.0
                + (x as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
            let new_y = model.plate_size.1 / 2.0
                + (y as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
            // println!("{}", (x + y) as f32%2.0);
            // l.push(Points::new(vec2(new_x, new_y), vec2(0.0, 0.0), vec2(0.0, 0.0), model.plate_size, (x + y) as f32%2.0 - 0.4));
            l.push(Points::new(
                vec2(new_x, new_y),
                vec2(0.0, 0.0),
                vec2(0.0, 0.0),
                model.plate_size,
                (x + y) as f32 % 2.0 - 0.5,
            ));
            // l.push(Points::new(vec2(new_x, new_y), vec2(0.0, 0.0), vec2(0.0, 0.0), model.plate_size, 0.5));
        }
    }
    model.p_l = l.clone();
    model.p_l_copy = l;
}

fn simulation_step(
    p_l_copy: &mut Vec<Points>,
    p_l: &mut Vec<Points>,
    app: &App,
    simulation_stop: &mut bool,
    model_total_force: &mut f32,
    epsilon: f32,
    sigma: f32,
) {
    *p_l_copy = p_l.clone();

    let mut total_force = 0.0;
    if *simulation_stop {
        return;
    }
    for (i, p) in p_l_copy.iter_mut().enumerate() {
        total_force += p.step(p_l, epsilon, sigma);
    }
    *model_total_force = total_force;

    *p_l = p_l_copy.clone();
}

fn main() {
    nannou::app(model).update(update).run();
}

fn model(app: &App) -> Model {
    let window_id = app
        .new_window()
        .view(view)
        .raw_event(raw_window_event)
        .build()
        .unwrap();
    let window = app.window(window_id).unwrap();
    let egui = Egui::from_window(&window);
    let num = 0.0;
    let seed_size = 100;
    let atom_count = 100;
    let mut model = Model {
        egui,
        sigma: 10.0,
        epsilon: 10.0,
        K: 20.0,
        plate_size: (300.0, 300.0),
        p_l_copy: vec![],
        p_l: vec![],
        speed: 1.0,
        time: 0.0,
        done_simulation: false,
        total_force: 0.0,
        seed_size,
        atom_count,
    };
    gen_atoms(&mut model, atom_count, seed_size);
    model
}

fn update(app: &App, model: &mut Model, update: Update) {
    let mut button_clicked = false;
    let mut slider_seed: Option<Response> = Option::None;
    let mut slider_atom: Option<Response> = Option::None;
    {
        let egui = &mut model.egui;
        egui.set_elapsed_time(update.since_start);

        let ctx = egui.begin_frame();

        egui::Window::new("Simulation controls").show(&ctx, |ui| {
            ui.label("step_size");
            ui.add(egui::Slider::new(&mut model.K, 1.0..=200.0));
            ui.label("sigma/distance");
            ui.add(egui::Slider::new(&mut model.sigma, 1.0..=100.0));
            ui.label("epsilon");
            ui.add(egui::Slider::new(&mut model.epsilon, 1.0..=1000.0));
            ui.label("speed");
            ui.add(egui::Slider::new(&mut model.speed, 1.0..=1000.0));
            ui.label("start control");
            ui.label("seed size");
            slider_seed = Some(ui.add(egui::Slider::new(&mut model.seed_size, 1..=500)));
            ui.label("atom count");
            slider_atom = Some(ui.add(egui::Slider::new(&mut model.atom_count, 1..=500)));

            ui.label("time");
            ui.label(model.time.to_string());
            ui.label("total potential");
            ui.label(model.total_force.to_string());

            let stop = ui.button("Stop");
            if stop.clicked() {
                model.done_simulation = !model.done_simulation;
                println!("stop: {}", model.done_simulation);
            }
            let step = ui.button("Step");
            if step.clicked() {
                for _ in 0..model.speed as u32 {
                    simulation_step(
                        &mut model.p_l_copy,
                        &mut model.p_l,
                        app,
                        &mut false,
                        &mut model.total_force,
                        model.epsilon,
                        model.sigma,
                    );
                    model.time += 0.01 * 0.1;
                }
            }

            let button = ui.button("Restart");
            button_clicked = button.clicked();
        });
    }
    if slider_seed.unwrap().changed() {
        gen_atoms(model, model.atom_count, model.seed_size);
    }
    if slider_atom.unwrap().changed() {
        gen_atoms(model, model.atom_count, model.seed_size);
    }
    if button_clicked {
        gen_atoms(model, model.atom_count, model.seed_size);
    }
    // if model.done_simulation {
    //     return;
    // }
    for _ in 0..model.speed as u32 {
        simulation_step(
            &mut model.p_l_copy,
            &mut model.p_l,
            app,
            &mut model.done_simulation,
            &mut model.total_force,
            model.epsilon,
            model.sigma,
        );
        if !model.done_simulation {
            model.time += 0.01 * 0.1;
            // break;
        }
    }
}

fn raw_window_event(_app: &App, model: &mut Model, event: &nannou::winit::event::WindowEvent) {
    model.egui.handle_raw_event(event);
}

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    draw.background().color(WHITE);
    // get app window size
    let size = app.window_rect();

    // draw circles
    for p in &model.p_l {
        // println!("p_l: {:?}", p);

        draw.ellipse()
            .x_y(
                p.pos.x - model.plate_size.0 / 2.0,
                p.pos.y - model.plate_size.1 / 2.0,
            )
            .radius(5.0)
            .color(Rgba::new(1.0, p.charge, 0.0, 1.0));
    }

    // draw.rect()
    //     .w_h(size.w(), size.h())
    //     .x_y(size.w() / 1.2, size.h() / 1.5)
    //     .stroke_color(BLACK)
    //     .stroke_weight(1.0)
    //     .no_fill();
    // draw_graph(&draw, app, &model);

    draw.to_frame(app, &frame).unwrap();
    model.egui.draw_to_frame(&frame).unwrap();
}
