use egui_plot::{Line, Plot, PlotPoints};
use nannou::{color::rgb::Rgb, draw::{properties::spatial::position, theme::Color}, prelude::*, text::rt::Point};
use nannou_egui::{
    self,
    egui::{self, lerp, Response},
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

 // TODO: Center the points around 0,0 insteada of 250 250
struct Model {
    // window: Window,
    egui: Egui,
    sigma: f32,             // nm
    epsilon: f32,           // J
    T: f32,                 // temperature
    plate_radius: f32,
    p_l: Vec<Points>,       // list of all particles
    p_l_copy: Vec<Points>,
    speed: f32,
    time: f32,
    stop_simulation: bool,
    avg_temperature: f32,
    seed_size: u32,
    atom_count: u32,
    temperature_depletion: f32,
}

fn gen_atoms(model: &mut Model, n: u32, seed_size: u32) {
    let mut l = vec![];
    
    let seed_x0 = (0.0 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
    let seed_y0 = (0.0 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
    let seed_x1 =
        (f32::sqrt(seed_size as f32) - 1.0) / 2.0 * 1.2 * model.sigma;
    let seed_y1 =
        (f32::sqrt(seed_size as f32) - 1.0) / 2.0 * 1.2 * model.sigma;
    for i in 0..n {
        let p = Points::random_init(model.plate_radius - 20.0, 0.5 - i as f32 % 2.0, model.T);
        if !(!(p.pos.x < seed_x0 - 20.0 || p.pos.x > seed_x1 + 20.0) && !(p.pos.y < seed_y0 - 20.0 || p.pos.y > seed_y1 + 20.0))
        {
            l.push(p);
        }
    }
    for x in 0..(f32::sqrt(seed_size as f32) as u32) {
        for y in 0..(f32::sqrt(seed_size as f32) as u32) {
            
            let new_x = (x as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
            let new_y = (y as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * 1.2 * model.sigma;
            // println!("{}", (x + y) as f32%2.0);
            // l.push(Points::new(vec2(new_x, new_y), vec2(0.0, 0.0), vec2(0.0, 0.0), model.plate_size, (x + y) as f32%2.0 - 0.4));
            l.push(Points::new(
                vec2(new_x, new_y),
                vec2(0.0, 0.0),
                vec2(0.0, 0.0),
                model.plate_radius,
                (x + y) as f32 % 2.0 - 0.5,
                true,
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
    model_avg_temperature: &mut f32,
    epsilon: f32,
    sigma: f32,
    temperature_depletion: f32,
) { // do a step of the simulation
    // *p_l_copy = p_l.clone(); // clone the list of particles so we can move everything at the same time

    let mut sum_temp = 0.0;
    if *simulation_stop { // to stop the simualtion
        return;
    }
    for p in p_l_copy.iter_mut() { // for each particle
        p.step(p_l, epsilon, sigma, temperature_depletion); // make a step
        sum_temp += p.temperature; // for the average temperature
    }
    *model_avg_temperature = sum_temp/p_l.len() as f32 * 1000.0; // calculate the average temperature

    *p_l = p_l_copy.clone(); // move the particles
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
    let seed_size = 100;
    let atom_count = 100;
    let mut model = Model {
        egui,
        sigma: 10.0,
        epsilon: 10.0,
        T: 20.0,
        temperature_depletion: 0.9995,
        plate_radius: 500.0,
        p_l_copy: vec![],
        p_l: vec![],
        speed: 1.0,
        time: 0.0,
        stop_simulation: false,
        avg_temperature: 0.0,
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
    let mut temp_slider: Option<Response> = Option::None;
    {
        let egui = &mut model.egui;
        egui.set_elapsed_time(update.since_start);

        let ctx = egui.begin_frame();

        egui::Window::new("Simulation controls").show(&ctx, |ui| {
            ui.label("temperature");
            temp_slider = Some(ui.add(egui::Slider::new(&mut model.T, 1.0..=100000.0)));
            ui.label("temperature depletion");
            ui.add(egui::Slider::new(&mut model.temperature_depletion, 0.999..=1.00000001));
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
            ui.label("average temperature");
            ui.label(model.avg_temperature.to_string());

            let stop = ui.button("Stop");
            if stop.clicked() { // stop simulation by toggling stop simulation boolean
                model.stop_simulation = !model.stop_simulation;
                println!("stop: {}", model.stop_simulation);
            }
            let step = ui.button("Step");
            if step.clicked() { // do a step of the simulation on click
                for _ in 0..model.speed as u32 {
                    simulation_step(
                        &mut model.p_l_copy,
                        &mut model.p_l,
                        app,
                        &mut false,
                        &mut model.avg_temperature,
                        model.epsilon,
                        model.sigma,
                        model.temperature_depletion,
                    );
                    model.time += 0.1;
                }
            }

            let button = ui.button("Restart");
            button_clicked = button.clicked();
        });
    }
    if slider_seed.unwrap().changed() || temp_slider.unwrap().changed() || slider_atom.unwrap().changed(){
        model.time = 0.0;
        gen_atoms(model, model.atom_count, model.seed_size);
    }
    if button_clicked {
        gen_atoms(model, model.atom_count, model.seed_size);
        model.time = 0.0;
    }
    for _ in 0..model.speed as u32 { // to speed up the simulation without sacrificing accuracy
        simulation_step(
            &mut model.p_l_copy,
            &mut model.p_l,
            app,
            &mut model.stop_simulation,
            &mut model.avg_temperature,
            model.epsilon,
            model.sigma,
            model.temperature_depletion,
        );
        if !model.stop_simulation {
            model.time += 0.1;
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
        if p.got_good{
            draw.ellipse()
                .x_y(
                    p.pos.x,
                    p.pos.y
                )
                .radius(5.0)
                .color(BLACK);
        } else {
            draw.ellipse()
                .x_y(
                    p.pos.x,
                    p.pos.y
                )
                .radius(5.0)
                .color(Hsl::new((p.charge+0.5)*250.0, 1.0, 0.5));

        }
    }
    draw.to_frame(app, &frame).unwrap();
    model.egui.draw_to_frame(&frame).unwrap();
}
