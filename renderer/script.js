var parsedData = [];
var canvas_width = 700;
var start_stop_button, reset_button, timestep_slider, timestep_label, timestep_value;
var num_bodies = 0;

function setup() {
  fetch('../data/output.csv')
    .then(response => response.text())
    .then(data => {
        // Parse the CSV data
        const rows = data.split('\n').slice(1);  // Split the CSV into rows and remove the header
        parsedData = rows.map(row => {
            if (row.trim() === '') return null;  // Skip blank lines
            const [time, mass, x, y] = row.split(',').map(Number);  // Split each row into columns and convert them to numbers
            return {time, mass, x, y};
        }).filter(Boolean);  // Remove null values from the array

        console.log(parsedData);  // Log the parsed data
        num_bodies = parsedData.filter(d => d.time == 0).length;
        print("Number of bodies: " + num_bodies);


        timestep_slider = select('#timestep-slider');
        timestep_label = select('#slider-label');
        timestep_value = timestep_slider.value();
        timestep_slider.attribute('max', parsedData[parsedData.length - 1].time);
        timestep_slider.input(() => {
          current_timestep = timestep_slider.value();
          timestep_slider.attribute('value', current_timestep);
          timestep_label.html(current_timestep);
          renderData();
          playing = false;
          start_stop_button.html('Start');
        });
        timestep_slider.attribute('value', current_timestep);
        timestep_label.html(current_timestep);

        let canvas = createCanvas(canvas_width, canvas_width);
        canvas.parent('canvas-container');
        background(20);
        stroke(255);
        noFill();
        frameRate(10);

        // Find the maximum position in the simulation in time 0
        //max_x_diff = max(parsedData.filter(d => d.time === 0).map(d => d.x)) - min(parsedData.filter(d => d.time === 0).map(d => d.x));
       // max_y_diff = max(parsedData.filter(d => d.time === 0).map(d => d.y)) - min(parsedData.filter(d => d.time === 0).map(d => d.y));
        //sim_pos_max = max(max_x_diff, max_y_diff);

        // Find maximum a*stddev(x) and a*stddev(y) for the simulation
        a=5;
        max_x_diff = a * math.std(parsedData.filter(d => d.time === 0).map(d => d.x));
        max_y_diff = a * math.std(parsedData.filter(d => d.time === 0).map(d => d.y));
        sim_pos_max = max(max_x_diff, max_y_diff);

        // Find the maximum mass in the simulation and calculate the mass scale
        max_mass = max(parsedData.filter(d => d.time === 0).map(d => d.mass));
        mass_scale = 16 / max_mass; // Scale for the mass of the bodies to be rendered, maximum mass is clamped at 16 solar masses

        // Scale the masses of the bodies
        for (let i = 0; i < parsedData.length; i++) {
          parsedData[i].mass *= mass_scale;
        }

        renderData();
    })
    .catch(error => console.error('Error:', error));

  

  start_stop_button = select('#start-stop');
  start_stop_button.mouseClicked(() => {
    playing = !playing;
    start_stop_button.html(playing ? 'Stop' : 'Start');
  });

  reset_button = select('#reset');
  reset_button.mouseClicked(() => {
    current_timestep = 0;
  });

  
}

var playing = false;
var current_timestep = 0;
var sim_pos_max; // Maximum position in the simulation (minimum is 0)
var sim_padding_scale = 0.1; // Padding around the simulation
const c = 2.99792e5; // Speed of light in km/s
const G = 4.3009e-3  // Gravitational constant in units of parsecs * (km/s)^2 / solar_mass
var mass_scale; // Scale for the mass of the bodies to be rendered
var trail = true; // Render the trail of the bodies

console.log("Canvas padding boundary: " + (canvas_width*sim_padding_scale) + " " + canvas_width*(1-sim_padding_scale))

function updateTimestepHTML() {
  timestep_slider.attribute('value', current_timestep);
  timestep_label.html(current_timestep);
}

function renderData() {
  updateTimestepHTML();

  background(20);

  // Slice parsedData where time is timestep
  
  let current_data = parsedData.filter(d => d.time === current_timestep);
  print("time: " + current_timestep + " data: " + current_data.length)
  
  for (let i = 0; i < current_data.length; i++) {
    //print(current_data[i].x + " " + current_data[i].y)
    let x = map(current_data[i].x, -0.5*sim_pos_max, 0.5*sim_pos_max, canvas_width*sim_padding_scale, canvas_width - canvas_width*sim_padding_scale);
    let y = map(current_data[i].y, -0.5*sim_pos_max, 0.5*sim_pos_max, canvas_width*sim_padding_scale, canvas_width - canvas_width*sim_padding_scale);
    let mass = current_data[i].mass;
    let r = map(mass, 0, 2, 1, 3, true);
    //fill(255, 255, 255, 80);
    if (mass <= 1.04) {
      colour = lerpColor(color(255, 153, 0, 80), color(255, 255, 0, 150), map(mass, 0.08, 1.04, 0, 1));
    } else if (mass <= 2.1) {
      colour = lerpColor(color(255, 255, 0, 150), color(255, 255, 255, 255), map(mass, 1.04, 2.1, 0, 1));
    } else {
      colour = lerpColor(color(255, 255, 255, 255), color(0, 255, 255, 255), map(mass, 2.1, 16, 0, 1));
    }
    fill(colour)
    noStroke();

    if (trail) { // If trail is enabled, render a line from the previous position to the current position
      let trail_data = parsedData.filter(d => d.time === current_timestep - 1);
      
      if (trail_data.length != 0) {
        let trail_x = map(trail_data[i].x, -0.5*sim_pos_max, 0.5*sim_pos_max, canvas_width*sim_padding_scale, canvas_width - canvas_width*sim_padding_scale);
        let trail_y = map(trail_data[i].y, -0.5*sim_pos_max, 0.5*sim_pos_max, canvas_width*sim_padding_scale, canvas_width - canvas_width*sim_padding_scale);

        stroke(colour);
        strokeWeight(2);
        line(trail_x, trail_y, x, y);
      } 
    }

    
    /*if (i == 0) {
      fill("orange")
      ellipse(x, y, 5);
      //r = (4 * G * mass) / c**2; // Schwarzschild radius
      r = 3;
      fill("black")
    }*/
    ellipse(x, y, r);
  }
}


function draw() {
  //print(playing)
  if (playing) {
    current_timestep++;
    if (current_timestep >= parsedData[parsedData.length - 1].time) {
      current_timestep = 0;
    }


    renderData()
  }
  
}