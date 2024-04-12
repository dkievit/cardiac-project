 /*-----------------------Start user-defined code ---------------------*/

 var sim; // Changing sim to var to enable "save image" button
 //model parameters, confusingly dubbed variables
 //for technical reasons
 var a = 0.025
 var e = 0.005
 var b = 4.5
 var dt = 3
 var disruption = 0.3

 let sumV; //to sum up V over the grid

 function cacatoo() {

   //setting details like simulation field size,
   //max run time, title, zoom scale, display
   //frequency, boundary conditions etc
   let config = {
     title: "Fitzhugh Nagumo PDE",
     description: "",
     maxtime: 1000000,
     ncol: 150,
     nrow: 150, // dimensions of the grid to build
     wrap: [false, false], // Wrap boundary [COLS, ROWS]   
     scale: 2, // scale of the grid (nxn pixels per grid cell
     bgcolour: 'white',
     graph_interval: 5,
     graph_update: 25
   }

   //make a model object
   sim = new Simulation(config)
   sim.makeGridmodel("Fitz");

   //set the colors for the 2D displays
   //arguments: which variable, max value, start color, end color
   sim.Fitz.colourGradient('V', 100, [255, 255, 255], [0, 0, 255])
   sim.Fitz.colourGradient('W', 100, [255, 255, 255], [255, 0, 0])

   //2D display for the V variable & colorbar
   sim.createDisplay_continuous({
     model: "Fitz",
     property: "V",
     label: "Action potential",
     minval: 0,
     maxval: 1,
     decimals: 3,
     nticks: 3
   })
   //2D display for the W variable & colorbar
   sim.createDisplay_continuous({
     model: "Fitz",
     property: "W",
     label: "Refractory tissue",
     minval: 0,
     maxval: 0.3,
     decimals: 3,
     nticks: 3
   })

  sim.reset = function(all = true) {
     if (all) sim.time = 0
     if (all) sim.Fitz.resetPlots()
     if (all) NewWaveAct = false
     if (all) AblateAct = false
     for (let i = 0; i < sim.Fitz.nc; i++) {
       for (let j = 0; j < sim.Fitz.nr; j++) {
         if (all) {
           if (i*i + j*j <= 100)
             sim.Fitz.grid[i][j].V = 1
           else
             sim.Fitz.grid[i][j].V = 0.0
           sim.Fitz.grid[i][j].W = 0.0

           sim.Fitz.grid[i][j].A = 0.025
         } else {
           if (i*i + j*j <= 100) sim.Fitz.grid[i][j].V = 1
         }
       }
     }
   }
   sim.reset()

	NewWaveAct = false
   sim.NewWave = function() {
     if (NewWaveAct) {
     for (let i = 0; i < sim.Fitz.nc; i++) {
       for (let j = 0; j < sim.Fitz.nr; j++) {
           if ((i-100)*(i-100) + (j-100)*(j-100) <= 100) sim.Fitz.grid[i][j].V = 1
         }
       }
     }
   }

	AblateAct = false
   sim.Ablate = function() {
     if (AblateAct) {
	     for (let i = 0; i < sim.Fitz.nc; i++) {
  	     for (let j = 0; j < sim.Fitz.nr; j++) {
    	       if ((i-100)*(i-100) + (j-100)*(j-100) <= 100) sim.Fitz.grid[i][j].V = 0
         }
       }
     }
   }

   sim.defibrillate = function() {
     for (let i = 0; i < sim.Fitz.nc; i++) {
       for (let j = 0; j < sim.Fitz.nr; j++) {
           sim.Fitz.grid[i][j].V = 1.0
           sim.Fitz.grid[i][j].W = 0.0   
       }
     }
   }
//   sim.defibrillate()


   //Defines the next-state function on a per position basis
   sim.Fitz.nextState = function(i, j) {
     //define variables V and W to store local V and W values
     let V = this.grid[i][j].V
     let W = this.grid[i][j].W
     let A = this.grid[i][j].A
     //update grid stored V and W values based on ODE equations
     this.grid[i][j].V += dt * (-V * (V - A) * (V - 1) - W)
     this.grid[i][j].W += dt * (e * (V - b * W))

   }

   //this is where the overall time loop of the simulation takes place
   sim.Fitz.update = function() {

     //this is for the simulation itself, resetting part of field
     if (sim.time % 100 == 0) sim.reset(false)
     if (sim.time % 50 == 0) sim.NewWave()
		 sim.Ablate()     
     
     //shifting of graphs, sliding window
     if (sim.time > 0 && sim.time % sim.config.graph_update * 10 == 0) {
       for (let name in sim.Fitz.graphs) {
         let g = sim.Fitz.graphs[name]
         let dat = g.data.slice(-sim.config.graph_update * 10)
         g.data = dat
         g.g.updateOptions({
           'file': this.data
         });
       }
     }

     //defines updating mode, calls update function on all positions in grid
     //causing the execution of only the ODE part
     this.synchronous()
     //now the PDE part of the equation is executed, note that only the V variable diffuses
     //second argument gives diffusion constant
     this.diffuseStates('V', 0.2)

     //here we compute the average value of V over the field
     sumV = 0
     for (let i = 0; i < sim.Fitz.nc; i++) {
       for (let j = 0; j < sim.Fitz.nr; j++) {
         sumV += this.grid[i][j].V
       }
     }
     
     //converting the sum to the average
     sumV /= sim.Fitz.nc * sim.Fitz.nr

     //this plots the V value in a particular position
     this.plotArray(["local V"],
       [this.grid[this.nc - 1][0].V],
       ["blue"],
       "Action potential at right-hand side")

     //this plots the field averaged V value
     this.plotArray(["average V"],
       [sumV],
       ["blue"],
       "Average action potential")

     //this adds an output screen showing the time and average V value
     //it is necessary to make this if you want to store this information to file (Save output data button)
     if (sim.time % 100 == 0) {
       if (sim.time == 0) sim.log("Time, avgV\n", "output")
       let line = `${sim.time}, `
       line += sumV + "\n"
       sim.log(line, "output")
     }
   } //end of time loop 

   //adding a button which when pressed makes a png picture of the 
   //V and W field at that time
   sim.addButton("Save grids", function() {
     sim.sectionToPNG("canvas_holder", "grid_timepoint_")
   })

   //adding a button which when pressed makes a png picture of the
   //two graphs (V in a point and avg V) we defined above at that time
   sim.addButton("Save graphs", function() {
     sim.sectionToPNG("graph_holder", "graph_timepoint_")
   })

   //adding a button which saves the output data in a flat text file
   //to enable plotting it in another program like e.g. gnuplot, excel
   sim.addButton("Save output data", function() {
     sim.write(document.getElementById('output').textContent, `Data_${sim.time}`)
   })

   sim.addHTML("form_holder", "<br>")

   sim.addButton("Refresh", function() {
     sim.reset()
   })

 	 sim.addButton("New Wave", function() {
     NewWaveAct = true
   })
   
   sim.addButton("Ablate", function() {
     AblateAct = true
   })

   sim.addButton("Defibrillate", function() {
     sim.defibrillate()
   })
   
   //adding the possibility of pulling a mouse over
   //the V screen puts V to zero there 
   sim.addStatebrush("Fitz", "V", 0.0, 40)

   //this starts the whole simulation, so it calls update
   sim.start()

 } //end overall cacatoo function


 /*-------------------------End user-defined code ---------------------*/
