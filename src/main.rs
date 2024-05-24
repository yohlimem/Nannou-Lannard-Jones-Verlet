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
}

fn gen_atoms(model: &mut Model, n: u32, seed_size: u32) {
    let mut l = vec![];
    // let points = model.plate_radius / 2.0;
    // make container
    // for i in 0..points as i32{
    //     let angle = lerp(0.0..=2.0*PI, i as f32/points as f32);
    //     let position = vec2(angle.cos(),angle.sin()) * model.plate_radius;
    //     l.push(Points::new(position, Vec2::ZERO, Vec2::ZERO, model.plate_radius, 0.0, false));
    // }

    for i in 0..n {
        l.push(Points::random_init(model.plate_radius - 20.0, 0.5 - i as f32 % 2.0, model.T))
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
) {
    *p_l_copy = p_l.clone();

    let mut sum_temp = 0.0;
    if *simulation_stop {
        return;
    }
    for p in p_l_copy.iter_mut() {
        p.step(p_l, epsilon, sigma);
        sum_temp += p.temperature;
    }
    *model_avg_temperature = sum_temp/p_l.len() as f32;

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
        T: 20.0,
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
            if stop.clicked() {
                model.stop_simulation = !model.stop_simulation;
                println!("stop: {}", model.stop_simulation);
            }
            let step = ui.button("Step");
            if step.clicked() {
                for _ in 0..model.speed as u32 {
                    simulation_step(
                        &mut model.p_l_copy,
                        &mut model.p_l,
                        app,
                        &mut false,
                        &mut model.avg_temperature,
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
    if slider_seed.unwrap().changed() || temp_slider.unwrap().changed() || slider_atom.unwrap().changed(){
        gen_atoms(model, model.atom_count, model.seed_size);
    }
    if button_clicked {
        gen_atoms(model, model.atom_count, model.seed_size);
    }
    for _ in 0..model.speed as u32 {
        simulation_step(
            &mut model.p_l_copy,
            &mut model.p_l,
            app,
            &mut model.stop_simulation,
            &mut model.avg_temperature,
            model.epsilon,
            model.sigma,
        );
        if !model.stop_simulation {
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
                p.pos.x,
                p.pos.y
            )
            .radius(5.0)
            .color(Hsl::new((p.charge+0.5)*250.0, 1.0, 0.5));
    }
    draw.to_frame(app, &frame).unwrap();
    model.egui.draw_to_frame(&frame).unwrap();
}
