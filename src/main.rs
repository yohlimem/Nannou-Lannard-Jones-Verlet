use egui_plot::{Line, Plot, PlotPoints};
use nannou::{color::rgb::Rgb, draw::theme::Color, prelude::*};
use nannou_egui::{self, egui, Egui};
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
}

fn gen_atoms(model: &mut Model, n: u32, seed_size: u32) {
    let mut l: Vec<Points> = Vec::new();
    let mut rng = rand::thread_rng();
    for i in 0..n {
        let x = rng.gen_range(0.0..model.plate_size.0);
        let y = rng.gen_range(0.0..model.plate_size.1);
        let p = Points::new(
            vec2(
                x,
                y,
            ),
            vec2(0.0, 0.0),
            vec2(0.0, 0.0),
            model.plate_size,
            i as f32%2.0  - 0.4,
        );
        // if p.pos.x
        //     < model.plate_size.0 / 2.0
        //         + ((f32::sqrt(seed_size as f32) as i32 - 1) / 2) as f32 * model.sigma * 2.0
        //         + 20.0
        //     && p.pos.x
        //         > model.plate_size.0 / 2.0
        //             - ((f32::sqrt(seed_size as f32) as i32) / 2) as f32 * model.sigma * 2.0
        //             - 20.0
        //     && p.pos.y
        //         < model.plate_size.1 / 2.0
        //             + ((f32::sqrt(seed_size as f32) as i32 - 1) / 2) as f32 * model.sigma * 2.0
        //             + 20.0
        //     && p.pos.y
        //         > model.plate_size.1 / 2.0
        //             - ((f32::sqrt(seed_size as f32) as i32) / 2) as f32 * model.sigma * 2.0
        //             - 20.0
        // {
        //     continue;
        // }
        l.push(p);
    }
    for x in 0..(f32::sqrt(seed_size as f32) as u32) {
        for y in 0..(f32::sqrt(seed_size as f32) as u32) {
            let new_x = model.plate_size.0 / 2.0 + (x as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * model.sigma * 0.3;
            let new_y = model.plate_size.1 / 2.0 + (y as f32 - ((seed_size as f32).sqrt() - 1.0) / 2.0) as f32 * model.sigma * 0.3;
            println!("{}", (x + y) as f32%2.0);
            // l.push(Points::new(vec2(new_x, new_y), vec2(0.0, 0.0), vec2(0.0, 0.0), model.plate_size, (x + y) as f32%2.0 - 0.4));
            l.push(Points::new(vec2(new_x, new_y), vec2(0.0, 0.0), vec2(0.0, 0.0), model.plate_size, (x + y) as f32%2.0 - 0.4));
        }
    }
    model.p_l = l.clone();
    model.p_l_copy = l;
}

fn simulation_step(p_l_copy: &mut Vec<Points>, p_l: &mut Vec<Points>) {
    // let mut new_p_l_copy = p_l.clone();
    *p_l_copy = p_l.clone();
    
    for p in p_l_copy.iter_mut() {
        p.step(p_l);
        // println!("step")
    }
    
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
    // let p_l: Vec<Points> = p_l.into_iter().map(|(x, y)| Points::new(vec2(x, y), vec2(0.0, 0.0), vec2(0.0, 0.0))).collect();
    let mut model = Model {
        egui,
        sigma: 100.0,
        epsilon: 500.0,
        K: 20.0,
        plate_size: (500.0, 500.0),
        p_l_copy: vec![],
        p_l: vec![],
        speed: 1.0,
    };
        gen_atoms(&mut model, 100, 50);
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
            ui.label("speed");
            ui.add(egui::Slider::new(&mut model.speed, 1.0..=100.0));
            // speed

            let button = ui.button("Restart");
            button_clicked = button.clicked();
        });
    }
    if button_clicked {
        gen_atoms(model, 100, 50);
    }
    for _ in 0..model.speed as u32 {
        simulation_step(&mut model.p_l_copy, &mut model.p_l);
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
        draw.ellipse()
            .x_y(
                p.pos.x - model.plate_size.0 / 2.0,
                p.pos.y - model.plate_size.1 / 2.0,
            )
            .radius(5.0)
            .color(Rgba::new(1.0, p.charge,0.0 , 1.0));
    }

    draw.rect()
        .w_h(size.w(), size.h())
        .x_y(size.w() / 1.2, size.h() / 1.5)
        .stroke_color(BLACK)
        .stroke_weight(1.0)
        .no_fill();
    // draw_graph(&draw, app, &model);

    draw.to_frame(app, &frame).unwrap();
    model.egui.draw_to_frame(&frame).unwrap();
}

// fn draw_graph(draw: &Draw, app: &App, model: &Model) {
//     let size = app.window_rect();
//     let graph_x = size.w() / 3.0;
//     let graph_y = size.h() / 3.0;
//     // let graph_x = 0.0;
//     // let graph_y = 0.0;

//     let size = vec2(3.0, 0.1);

//     for r in -5..1000 {
//         if r == 0 || r == 1 {
//             continue;
//         }
//         let start_x = (r as f32 - 1.0) * size.x + graph_x;
//         let start_y = u(&model, r as f32 - 1.0) * size.y + graph_y;
//         let end_x = r as f32 * size.x + graph_x;
//         let end_y = u(&model, r as f32) * size.y + graph_y;

//         draw.line()
//             .start(pt2(start_x, start_y))
//             .end(pt2(end_x, end_y))
//             .color(BLACK);
//     }
// }
