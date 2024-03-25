use nannou::prelude::*;
use nannou_egui::{self, egui, Egui};
use rand::Rng;
use egui_plot::{Line, Plot, PlotPoints};

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
    sigma: f32, // nm
    epsilon: f32, // J
    K: f32, // temporary
    plate_size: (f32, f32), // nm^2
    p_l: Vec<(f32, f32)>, // list of all particles
    p_l_copy: Vec<(f32, f32)>,
}

fn gen_atoms(model: &mut Model, n: u32, seed_size: u32){
    let mut l:Vec<(f32, f32)> = Vec::new();
    let mut rng = rand::thread_rng();
    for _ in 0..n{
        let p = (rng.gen_range(0.0..model.plate_size.0 as f32), rng.gen_range(0.0..model.plate_size.1 as f32));
        if p.0 < model.plate_size.0/2.0+((f32::sqrt(seed_size as f32) as i32-1)/2) as f32 * model.sigma*2.0 + 20.0 && p.0 > model.plate_size.0/2.0 - ((f32::sqrt(seed_size as f32) as i32)/2) as f32 * model.sigma*2.0 - 20.0 && p.1 < model.plate_size.1/2.0+((f32::sqrt(seed_size as f32) as i32-1)/2) as f32 * model.sigma*2.0 + 20.0 && p.1 > model.plate_size.1/2.0 - ((f32::sqrt(seed_size as f32) as i32)/2) as f32 * model.sigma*2.0 - 20.0 {
            continue;
        }
        l.push(p);
    }
    for x in 0..(f32::sqrt(seed_size as f32) as u32){
        for y in 0..(f32::sqrt(seed_size as f32) as u32){
            let new_x = model.plate_size.0 / 2.0 + (x as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * model.sigma * 2.0;
            let new_y = model.plate_size.1 / 2.0 + (y as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * model.sigma * 2.0;

            l.push((new_x, new_y))
        }
    }
    model.p_l = l.clone();
    model.p_l_copy = l;
    // return l;
    // p_l = gen_atoms(200, 50);
}
fn u(model: &Model, r: f32) -> f32 { // potential energy
    return 4.0 * model.epsilon * ((model.sigma / r).powi(12) - (model.sigma / r).powi(6));
}

fn distance(model: &Model, i: usize, j: usize) -> f32 {
    let p_i = model.p_l[i];
    let p_j = model.p_l[j];
    ((p_j.0 - p_i.0).powi(2) + (p_j.1 - p_i.1).powi(2)).sqrt()
}

fn distance_pos(p: (f32, f32), p2: (f32, f32)) -> f32 {
    ((p2.0 - p.0).powi(2) + (p2.1 - p.1).powi(2)).sqrt()
}

fn energy(model: &Model, i: usize) -> f32 {
    let mut potential = 0.0;
    for j in 0..model.p_l.len() {
        if j == i {
            continue;
        }
        potential += u(&model, distance(model, i, j));
    }
    potential
}

fn energy_pos(model: &Model, p: (f32, f32), index: usize) -> f32 {
    let mut potential = 0.0;
    for (i, p2) in model.p_l.iter().enumerate() {
        if i == index {
            continue;
        }
        potential += u(&model, distance_pos(p, *p2));
    }
    potential
}

fn choose_delta(model: &Model, i: usize) -> Option<(f32, f32)> {
    let pos = model.p_l[i];
    for _ in 0..10 {
        let dy = rand::thread_rng().gen_range(-model.K..model.K);
        let dx = rand::thread_rng().gen_range(-model.K..model.K);
        if pos.0 + dx < model.plate_size.0 && pos.0 + dx > 0.0 && pos.1 + dy > 0.0 || pos.1 + dx < model.plate_size.1 {
            return Some((pos.0 + dx, pos.1 + dy));
        }
    }
    None
}

fn step(model: &mut Model, i: usize) {
    if let Some(new_pos) = choose_delta(model, i) {
        if energy_pos(model, new_pos, i) < energy(model, i) {
            model.p_l_copy[i] = new_pos;
        }
    }
}

fn simulation(model: &mut Model) {
    for i in 0..model.p_l.len() {
        step(model, i);
    }
    model.p_l = model.p_l_copy.clone();
}

// fn simulation_loop(model: &mut Model, n_iter: usize) {
//     for _ in 0..n_iter {
//         simulation(model);
//     }
// }

fn main() {
    nannou::app(model).update(update).run();
    
}

fn model(app: &App) -> Model {
    let window_id = app.new_window().view(view).raw_event(raw_window_event).build().unwrap();
    let window = app.window(window_id).unwrap();
    let egui = Egui::from_window(&window);
    let num = 0.0;
    let p_l = vec![(100.0, 100.0), (60.0, 35.0), (20.0, 200.0), (300.0, 300.0)];
    let mut model = Model {
        egui,
        sigma: 10.0,
        epsilon: 500.0,
        K: 20.0,
        plate_size: (500.0, 500.0),
        p_l_copy: p_l.clone(),
        p_l,
    };
    gen_atoms(&mut model, 100, 20);
    // println!("{:?}", model.p_l);
    model
}

fn update(app: &App, model: &mut Model, update: Update) {
    let mut button_clicked = false;
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
            ui.label("Restart simulation");

            
            let button = ui.button("Restart");
            button_clicked = button.clicked()
        });

    }
    if button_clicked {
        gen_atoms(model, 100, 20);
    }

    simulation(&mut *model);


}

fn raw_window_event(_app: &App, model: &mut Model, event: &nannou::winit::event::WindowEvent){
    model.egui.handle_raw_event(event);
}

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    draw.background().color(WHITE);
    // get app window size
    let size = app.window_rect();


    // draw circles
    for p in &model.p_l {
        draw.ellipse().x_y(p.0 - model.plate_size.0 / 2.0, p.1 - model.plate_size.1 / 2.0).radius(5.0).color(BLACK);
    }

    draw.rect().w_h(size.w(), size.h()).x_y(size.w()/1.2, size.h()/1.5).stroke_color(BLACK).stroke_weight(1.0).no_fill();
    draw_graph(&draw, app, &model);

    
    draw.to_frame(app, &frame).unwrap();
    model.egui.draw_to_frame(&frame).unwrap();
}

fn draw_graph(draw: &Draw, app: &App, model: &Model){
    let size = app.window_rect();
    let graph_x = size.w() / 3.0;
    let graph_y = size.h() / 3.0;
    // let graph_x = 0.0;
    // let graph_y = 0.0;
    
    let size = vec2(3.0, 0.1);

    for r in -5..1000{
        if r == 0 || r == 1{
            continue;
        }
        let start_x = (r as f32 - 1.0) * size.x + graph_x;
        let start_y = u(&model, r as f32 - 1.0) * size.y + graph_y;
        let end_x = r as f32 * size.x + graph_x;
        let end_y = u(&model, r as f32) * size.y + graph_y;

        draw.line()
            .start(pt2(start_x, start_y))
            .end(pt2(end_x, end_y))
            .color(BLACK);
    }
}
