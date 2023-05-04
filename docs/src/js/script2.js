function colors(){ 
    var red= document.getElementById("red").value; 
    var green = document.getElementById("green").value; 
    var blue = document.getElementById("blue").value; 
    return red, green, blue;
}

console.log("boombeach")

/* When the user clicks on the button, 
toggle between hiding and showing the dropdown content */
function myFunction() {
    document.getElementById("myDropdown").classList.toggle("show");
} 
function myFunction2() {
    document.getElementById("myDropdown2").classList.toggle("show");
}
function myFunction3() {
    document.getElementById("myDropdown3").classList.toggle("show");
}

// Close the dropdown if the user clicks outside of it
window.onclick = function(event) {
    if (!event.target.matches('.dropbtn')) {
        var dropdowns = document.getElementsByClassName("dropdown-content");
        var i;
        for (i = 0; i < dropdowns.length; i++) {
            var openDropdown = dropdowns[i];
            if (openDropdown.classList.contains('show')) {
                openDropdown.classList.remove('show');
            }
        }
    }
}

var canvas = document.getElementById("myCanvas");
var c = canvas.getContext("2d");
canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 100;

canvas.focus();	

var simHeight = 1.14;	
var cScale = canvas.height / 1.0;
var simWidth = canvas.width / cScale;

var U_FIELD = 0;
var V_FIELD = 1;
var S_FIELD = 2;

var cnt = 0;

function cX(x) {
    return x * cScale;
}

function cY(y) {	
    return canvas.height - y * cScale;
}

// ----------------- start of simulator ------------------------------

class Fluid1 {
    constructor(density, numX, numY, h) {
        this.density = density;
        this.numX = numX + 2; 
        this.numY = numY + 2;
        this.numCells = this.numX * this.numY;
        this.h = h;
        this.u = new Float32Array(this.numCells);
        this.u1 = new Float32Array(this.numCells);
        this.v = new Float32Array(this.numCells);
        this.v1 = new Float32Array(this.numCells);
        this.newU = new Float32Array(this.numCells);
        this.newU1 = new Float32Array(this.numCells);
        this.newV1 = new Float32Array(this.numCells);
        this.newV = new Float32Array(this.numCells);
        this.p = new Float32Array(this.numCells);
        this.s = new Float32Array(this.numCells);
        this.m = new Float32Array(this.numCells);
        this.m1 = new Float32Array(this.numCells);
        this.newM = new Float32Array(this.numCells);
        this.newM1 = new Float32Array(this.numCells);
        this.m.fill(1.0);
        this.m1.fill(1.0);
        var num = numX * numY;
    }

    integrate(dt, gravity) {
        var n = this.numY;
        for (var i = 1; i < this.numX; i++) {
            for (var j = 1; j < this.numY-1; j++) {
                if (this.s[i*n + j] != 0.0 && this.s[i*n + j-1] != 0.0) {
                    this.v[i*n + j] += gravity * dt;
                }
            }	 
        }
    }

    solveIncompressibility(numIters, dt) {

        var n = this.numY;
        var cp = this.density * this.h / dt;

        for (var iter = 0; iter < numIters; iter++) {

            for (var i = 1; i < this.numX-1; i++) {
                for (var j = 1; j < this.numY-1; j++) {
                    // if (this.s[i*n + j] == 0.0)
                    // 	continue;
                    if (this.s[i*n + j] != 0.0) {
                        var s = this.s[i*n + j];
                        var sx0 = this.s[(i-1)*n + j];
                        var sx1 = this.s[(i+1)*n + j];
                        var sy0 = this.s[i*n + j-1];
                        var sy1 = this.s[i*n + j+1];
                        var s = sx0 + sx1 + sy0 + sy1;
                        // if (s == 0.0)
                        // 	continue;
                        if (s != 0.0) {
                            var div = this.u[(i+1)*n + j] - this.u[i*n + j] + 
                                this.v[i*n + j+1] - this.v[i*n + j];
                            var div1 = this.u1[(i+1)*n + j] - this.u1[i*n + j] + 
                            this.v1[i*n + j+1] - this.v1[i*n + j];
                            var p = -div / s;
                            p *= scene.overRelaxation;
                            this.p[i*n + j] += cp * p;
                            var p1 = -div1/s;
                            p1 *= scene.overRelaxation;
                            this.u[i*n + j] -= sx0 * p;
                            this.u[(i+1)*n + j] += sx1 * p;
                            this.v[i*n + j] -= sy0 * p;
                            this.v[i*n + j+1] += sy1 * p;
                            this.u1[i*n + j] -= sx0 * p1;
                            this.u1[(i+1)*n + j] += sx1 * p1;
                            this.v1[i*n + j] -= sy0 * p1;
                            this.v1[i*n + j+1] += sy1 * p1;
                        }
                    }
                        
                }
            }
        }
    }

    grenadeCompressibility(numIters, dt) {

        var n = this.numY;
        var cp = this.density * this.h / dt;
        // console.log("True");
            for (var i = 1; i < this.numX-1; i++) {
                for (var j = 1; j < this.numY-1; j++) {
                    // if (this.s[i*n + j] == 0.0)
                    // 	continue;
                    if(this.grenadeX == i && this.grenadeY == j){
                        this.p[i*n + j] = 10;
                    }
            }
        }
    }

    extrapolate() {
        var n = this.numY;
        for (var i = 0; i < this.numX; i++) {
            this.u[i*n + 0] = this.u[i*n + 1];
            this.u[i*n + this.numY-1] = this.u[i*n + this.numY-2]; 
            this.u1[i*n + 0] = this.u1[i*n + 1];
            this.u1[i*n + this.numY-1] = this.u1[i*n + this.numY-2]; 
        }
        for (var j = 0; j < this.numY; j++) {
            this.v[0*n + j] = this.v[1*n + j];
            this.v[(this.numX-1)*n + j] = this.v[(this.numX-2)*n + j];
            this.v1[0*n + j] = this.v1[1*n + j];
            this.v1[(this.numX-1)*n + j] = this.v1[(this.numX-2)*n + j];
        }
    }

    sampleField(x, y, field) {
        var n = this.numY;
        var h = this.h;
        var h1 = 1.0 / h;
        var h2 = 0.5 * h;

        x = Math.max(Math.min(x, this.numX * h), h);
        y = Math.max(Math.min(y, this.numY * h), h);

        var dx = 0.0;
        var dy = 0.0;

        var f;
        var f1;

        switch (field) {
            case U_FIELD: f = this.u; f1 = this.u1; dy = h2; break;
            case V_FIELD: f = this.v; f1 = this.v1; dx = h2; break;
            case S_FIELD: f = this.m; dx = h2; dy = h2; f1 = this.m1; break

        }

        var x0 = Math.min(Math.floor((x-dx)*h1), this.numX-1);
        var tx = ((x-dx) - x0*h) * h1;
        var x1 = Math.min(x0 + 1, this.numX-1);
        
        var y0 = Math.min(Math.floor((y-dy)*h1), this.numY-1);
        var ty = ((y-dy) - y0*h) * h1;
        var y1 = Math.min(y0 + 1, this.numY-1);

        var sx = 1.0 - tx;
        var sy = 1.0 - ty;

        var val = sx*sy * f[x0*n + y0] +
            tx*sy * f[x1*n + y0] +
            tx*ty * f[x1*n + y1] +
            sx*ty * f[x0*n + y1];
        var val1 = sx*sy * f1[x0*n + y0] +
            tx*sy * f1[x1*n + y0] +
            tx*ty * f1[x1*n + y1] +
            sx*ty * f1[x0*n + y1];
        return [val, val1];
    }

    avgU(i, j) {
        var n = this.numY;
        var u = (this.u[i*n + j-1] + this.u[i*n + j] +
            this.u[(i+1)*n + j-1] + this.u[(i+1)*n + j]) * 0.25;
        return u;
    }

    avgU1(i, j) {
        var n = this.numY;
        var u = (this.u1[i*n + j-1] + this.u1[i*n + j] +
            this.u1[(i+1)*n + j-1] + this.u1[(i+1)*n + j]) * 0.25;
        return u;
    }

    avgV(i, j) {
        var n = this.numY;
        var v = (this.v[(i-1)*n + j] + this.v[i*n + j] +
            this.v[(i-1)*n + j+1] + this.v[i*n + j+1]) * 0.25;
        return v;
    }

    avgV1(i, j) {
        var n = this.numY;
        var v = (this.v1[(i-1)*n + j] + this.v1[i*n + j] +
            this.v1[(i-1)*n + j+1] + this.v1[i*n + j+1]) * 0.25;
        return v;
    }


    advectVel(dt) {

        this.newU.set(this.u);
        this.newV.set(this.v);
        this.newU1.set(this.u1);
        this.newV1.set(this.v1);

        var n = this.numY;
        var h = this.h;
        var h2 = 0.5 * h;

        for (var i = 1; i < this.numX; i++) {
            for (var j = 1; j < this.numY; j++) {

                cnt++;

                // u component
                if (this.s[i*n + j] != 0.0 && this.s[(i-1)*n + j] != 0.0 && j < this.numY - 1) {
                    var x = i*h;
                    var y = j*h + h2;
                    var x1 = x;
                    var y1 = y;
                    var u = this.u[i*n + j];
                    var u1 = this.u1[i*n + j];
                    var v = this.avgV(i, j);
                    var v1 = this.avgV1(i, j);
//						var v = this.sampleField(x,y, V_FIELD);
                    x = x - dt*u;
                    x1 = x1 - dt*u1;
                    y = y - dt*v;
                    y1 = y1 - dt*v1;
                    u = this.sampleField(x,y, U_FIELD);
                    u1 = this.sampleField(x1, y1, U_FIELD);
                    this.newU[i*n + j] = u[0];
                    this.newU1[i*n + j] = u1[1];
                }
                // v component
                if (this.s[i*n + j] != 0.0 && this.s[i*n + j-1] != 0.0 && i < this.numX - 1) {
                    var x = i*h + h2;
                    var y = j*h;
                    var x1 = x;
                    var y1 = y;
                    var u = this.avgU(i, j);
                    var u1 = this.avgU1(i, j);
//						var u = this.sampleField(x,y, U_FIELD);
                    var v = this.v[i*n + j];
                    var v1 = this.v1[i*n + j];
                    x = x - dt*u;
                    x1 = x1 - dt*u1;
                    y = y - dt*v;
                    y1 = y1 - dt*v1;
                    v = this.sampleField(x,y, V_FIELD);
                    v1 = this.sampleField(x1, y1, V_FIELD);
                    this.newV[i*n + j] = v[0];
                    this.newV1[i*n + j] = v1[1];
                }
            }	 
        }

        this.u.set(this.newU);
        this.v.set(this.newV);
        this.u1.set(this.newU1);
        this.v1.set(this.newV1);
    }

    advectSmoke(dt) {

        this.newM.set(this.m);
        this.newM1.set(this.m1);	
        var n = this.numY;
        var h = this.h;
        var h2 = 0.5 * h;

        for (var i = 1; i < this.numX-1; i++) {
            for (var j = 1; j < this.numY-1; j++) {

                if (this.s[i*n + j] != 0.0) {
                    var u = (this.u[i*n + j] + this.u[(i+1)*n + j]) * 0.5;
                    var v = (this.v[i*n + j] + this.v[i*n + j+1]) * 0.5;
                    var x = i*h + h2 - dt*u;
                    var y = j*h + h2 - dt*v;

                    this.newM[i*n + j] = this.sampleField(x,y, S_FIELD)[0];
                    var u1 = (this.u1[i*n + j] + this.u1[(i+1)*n + j]) * 0.5;
                    var v1 = (this.v1[i*n + j] + this.v1[i*n + j+1]) * 0.5;
                    var x1 = i*h + h2 - dt*u1;
                    var y1 = j*h + h2 - dt*v1;
                    this.newM1[i*n + j] = this.sampleField(x1, y1, S_FIELD)[1];
                 }
            }	 
        }
        this.m.set(this.newM);
        this.m1.set(this.newM1);
    }

    curl(x, y) {
        var n = this.numY;
        // console.log(this);
        return this.u[x * n + y + 1] - this.u[x * n + y - 1] + this.v[(x - 1) * n + y] + this.v[(x + 1) * n + y];
    }

    vorticityConfinement(dt) {

        this.newU.set(this.u);
        this.newV.set(this.v);

        var width = this.numX;
        var height = this.numY;

        var dx = 0.0;
        var dy = 0.0;
        var len = 0.0;
        var vorticity = 10.0; //10.0

        var n = this.numY;
        var h = this.h;
        var h2 = 0.5 * h;

        for (var my_y = 2; my_y < this.numY - 3; my_y++) {
            for (var my_x = 2; my_x < this.numX - 3; my_x++) {
                dx = Math.abs(this.curl(my_x + 0, my_y - 1)) - Math.abs(this.curl(my_x + 0, my_y + 1));
                dy = Math.abs(this.curl(my_x + 1, my_y + 0)) - Math.abs(this.curl(my_x - 1, my_y + 0));
                len = Math.sqrt(dx ** 2 + dy ** 2) + 10 ** -5;
                dx = vorticity / len * dx;
                dy = vorticity / len * dy;
                this.newU[my_x * n + my_y] = this.u[my_x * n + my_y] + dt * this.curl(my_x, my_y) * dx;
                this.newV[my_x * n + my_y] = this.v[my_x * n + my_y] + dt * this.curl(my_x, my_y) * dy;
            }
        }
        this.u.set(this.newU);
        this.v.set(this.newV);
    }


    // ----------------- end of simulator ------------------------------


    simulate(dt, gravity, numIters) {

        if(this.showGrenade){
            this.integrate(dt, gravity);
            this.p.fill(0,0);
            this.grenadeCompressibility(x, y);
            this.extrapolate();
            this.advectVel(dt);
            this.vorticityConfinement(dt);
            this.advectSmoke(dt);
        }

        this.integrate(dt, gravity);

        this.p.fill(0.0);
        this.solveIncompressibility(numIters, dt);

        this.extrapolate();
        this.advectVel(dt);
        this.vorticityConfinement(dt);
        this.advectSmoke(dt);
    }
}

var scene = 
{
    gravity : -9.81,
    dt : 1.0 / 120.0,
    numIters : 100,
    frameNr : 0,
    overRelaxation : 1.9,
    obstacleX : 0.0,
    obstacleY : 0.0,
    obstacleRadius: 0.15,
    paused: false,
    sceneNr: 0,
    showObstacle: false,
    showStreamlines: false,
    showVelocities: false,	
    showPressure: false,
    showSmoke: true,
    showGrenade: false,
    realistic_shader: false,
    grenadeX: 0.0,
    grenadeY: 0.0,
    fluid: null, 
    time: 0.0,
    activeSmokes: [], 
    activeFlares: [],
    smokegrenadeRadius: 0.1,
    flaregrenadeRadius: 0.01,
    smokegrenadeDuration: 3.0, 
    flaregrenadeDuration: 1.0, 
};

function checkScene(){
    if(scene.boombeach){
        loadScript("js/script.js");
    }
    
    function loadScript(src){
        var el = document.createElement("script");
        el.src = src;
        document.body.appendChild(el);
    }
}

function setupScene(sceneNr = 0) 
{
    scene.sceneNr = sceneNr;
    scene.obstacleRadius = 0.15;
    scene.grenadeRadius = 0.1;
    scene.overRelaxation = 1.9;

    scene.dt = 1.0 / 60.0;
    scene.numIters = 40;

    var res = 100;
    
    if (sceneNr == 0)
        res = 50;
    else if (sceneNr == 3)
        res = 200;

    var domainHeight = 1.0;
    var domainWidth = domainHeight / simHeight * simWidth;
    var h = domainHeight / res;

    var numX = Math.floor(domainWidth / h);
    var numY = Math.floor(domainHeight / h);

    var density = 1000.0;

    f = scene.fluid = new Fluid1(density, numX, numY, h);

    var n = f.numY;

    if (sceneNr == 0) {   		// tank

        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {
                var s = 1.0;	// fluid
                if (i == 0 || i == f.numX-1 || j == 0)
                    s = 0.0;	// solid
                f.s[i*n + j] = s
            }
        }
        scene.gravity = -9.81;
        scene.showPressure = true;
        scene.showSmoke = false;
        scene.showStreamlines = false;
        scene.showVelocities = false;
        scene.ball = true; //added
        scene.explosion = false; //added
        scene.flaregrenade = false; //added
        scene.smokegrenade = false; //added
        scene.showTemp = false; // added
        scene.smokeAppOne = true; // added
        scene.smokeAppTwo = false; // added
        scene.smokeAppThree = false; // added
        scene.fire = false;
    }
    else if (sceneNr == 1 || sceneNr == 3) { // vortex shedding

        var inVel = 2.0;
        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {
                var s = 1.0;	// fluid
                if (i == 0 || j == 0 || j == f.numY-1)
                    s = 0.0;	// solid
                f.s[i*n + j] = s

                if (i == 1) {
                    f.u[i*n + j] = inVel;
                }
            }
        }

        var pipeH = 0.1 * f.numY;
        var minJ = Math.floor(0.5 * f.numY - 0.5*pipeH);
        var maxJ = Math.floor(0.5 * f.numY + 0.5*pipeH);

        for (var j = minJ; j < maxJ; j++){
            f.m[j] = 0.0;
            f.m1[j] = 0.0;
        }


        if (scene.ball) {
            setObstacle(0.4, 0.5, true)
        }
    
        scene.gravity = 0.0;
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;
        scene.ball = false; 
        scene.showTemp = false;
        scene.smokeAppOne = true; // added
        scene.smokeAppTwo = false; // added
        scene.smokeAppThree = false; // added
        scene.fire = false;

        if (sceneNr == 3) {
            scene.dt = 1.0 / 120.0;
            scene.numIters = 100;
            scene.showPressure = true;
        }

    }
    else if (sceneNr == 2) { // paint

        scene.gravity = 0.0;
        scene.overRelaxation = 1.0;
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;
        scene.obstacleRadius = 0.1;
        scene.ball = false; 
        scene.showTemp = false;
    scene.smokeAppOne = true; // added
    scene.smokeAppTwo = false; // added
    scene.smokeAppThree = false; // added
    scene.fire = false;
    } else if (sceneNr == 4) { // None

        scene.gravity = 0.0;
        scene.overRelaxation = 1.0;
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;
        scene.obstacleRadius = 0.1;
        scene.ball = false; 
        scene.activeFlares = []
        scene.showTemp = false;
    scene.smokeAppOne = true; // added
    scene.smokeAppTwo = false; // added
    scene.smokeAppThree = false; // added
    scene.fire = false;
    }

    document.getElementById("smokeButton").checked = scene.showSmoke;
    document.getElementById("overrelaxButton").checked = scene.overRelaxation > 1.0;
    document.getElementById("ballButton").checked = scene.ball;
    document.getElementById("explosionButton").checked = scene.explosion;
    document.getElementById("flaregrenadeButton").checked = scene.flaregrenade;
    document.getElementById("smokegrenadeButton").checked = scene.smokegrenade;
    document.getElementById("tempButton").checked = scene.showTemp;
    document.getElementById("smokeAppOne").checked = scene.smokeAppOne;
    document.getElementById("smokeAppTwo").checked = scene.smokeAppTwo;
    document.getElementById("smokeAppThree").checked = scene.smokeAppThree;
    document.getElementById("fire").checked = scene.fire;
}

// draw -------------------------------------------------------

function setColor(r,g,b) {
    c.fillStyle = `rgb(
        ${Math.floor(255*r)},
        ${Math.floor(255*g)},
        ${Math.floor(255*b)})`
    c.strokeStyle = `rgb(
        ${Math.floor(255*r)},
        ${Math.floor(255*g)},
        ${Math.floor(255*b)})`
}

function getSciColor(val, minVal, maxVal) {
    val = Math.min(Math.max(val, minVal), maxVal- 0.0001);
    var d = maxVal - minVal;
    val = d == 0.0 ? 0.5 : (val - minVal) / d;
    var m = 0.25;
    var num = Math.floor(val / m);
    var s = (val - num * m) / m;
    var r, g, b;

    switch (num) {
        case 0 : r = 0.0; g = s; b = 1.0; break;
        case 1 : r = 0.0; g = 1.0; b = 1.0-s; break;
        case 2 : r = s; g = 1.0; b = 0.0; break;
        case 3 : r = 1.0; g = 1.0 - s; b = 0.0; break;
    }

    return[255*r,255*g,255*b, 0]
}

function norm(vector) {
    var vectorLength = vector.length;
    var total = 0.0;
    for (var i = 0; i < vectorLength; i++) {
        total += vector[i] * vector[i];
    }
    // console.log("total: ", total)
    return Math.sqrt(total);
}

function normalize(vector) {
    var vectorNorm = norm(vector);
    if (vectorNorm == 0) {
        return vector;
    }
    var finalVector = []
    var vectorLength = vector.length;
    for (var i = 0; i < vectorLength; i++) {
        finalVector.push(vector[i] / vectorNorm);
    }
    // console.log("finalVector: ", finalVector);
    return finalVector;
}

function dot(vector1, vector2) {
    var vectorLength = vector1.length;
    total = 0.0
    for (var i = 0; i < vectorLength; i++) {
        total += vector1[i] * vector2[i];
    }
    return total;
}

function lightingShader(mass, velocity, ambient=0.2) {
    var normVelocity = normalize(velocity);
    var normLightDir = normalize([1.0, 1.0]);
    var store = dot(normVelocity, normLightDir);
    // console.log("store: ", store);
    
    var shading = Math.max(store, 0.0) + ambient;
    shading = Math.min(shading, 1.0);
    return shading
}

function realisticShader(id, f){
    var cx = Math.floor(cScale * 1.1 * f.h) + 1;
    var cy = Math.floor(cScale * 1.1 * f.h) + 1;
    var copy = id;
    for (var i = 0; i < f.numX; i++){
        for (var j = 0; j < f.numY; j++){
            var x = Math.floor(cX(i * f.h));
            var y = Math.floor(cY((j+1) * f.h));
            for (var yi = y; yi < y + cy; yi++) {
                var p = 4 * (yi * canvas.width + x);
                for (var xi = 0; xi < cx; xi++) {
                    if(p - 3 >= 0 && id.length > p+4*canvas.width +7 && p - 4*canvas.width - 3 >= 0){
                        avg_r = id[p+1] + id[p-3] + id[p+5] + id[p+4*(canvas.width) + 1] + id[p - 4*canvas.width + 1] + id[p+4*(canvas.width) + 5] + id[p+4*(canvas.width) - 3] + id[p - 4*canvas.width + 5] + + id[p - 4*canvas.width - 3];
                        avg_g = id[p+2] + id[p-2] + id[p+6] + id[p+4*(canvas.width) + 2] + id[p - 4*canvas.width + 2] + id[p+4*(canvas.width) + 6] + id[p+4*(canvas.width) - 2] + id[p - 4*canvas.width + 6] + + id[p - 4*canvas.width - 2];
                        avg_b = id[p+3] + id[p-1] + id[p+7] + id[p+4*(canvas.width) + 3] + id[p - 4*canvas.width + 3] + id[p+4*(canvas.width) + 7] + id[p+4*(canvas.width) - 1] + id[p - 4*canvas.width + 7] + + id[p - 4*canvas.width - 1];
                        avg_a = 
                        copy[p++] = avg_r/9;
                        copy[p++] = avg_g/9;
                        copy[p++] = avg_b/9;
                        copy[p++] = id[p+4];
                    } else {
                        copy[p++] = id[p++];
                        copy[p++] = id[p++];
                        copy[p++] = id[p++];
                        copy[p++] = id[p++];
                    }
                }
            }
        }
    }
    id = copy;
    return id;
}

function draw() 
{
    c.clearRect(0, 0, canvas.width, canvas.height);

    c.fillStyle = "#FF0000";
    f = scene.fluid;
    n = f.numY;

    var cellScale = 1.1;

    var h = f.h;

    minP = f.p[0];
    maxP = f.p[0];
    for (var i = 0; i < f.numCells; i++) {
        minP = Math.min(minP, f.p[i]);
        maxP = Math.max(maxP, f.p[i]);
    }

    id = c.getImageData(0,0, canvas.width, canvas.height)

    var color = [255, 255, 255, 100]
    var color1 = [255, 255, 255, 100]

    for (var i = 0; i < f.numX; i++) {
        for (var j = 0; j < f.numY; j++) {

            if (scene.showPressure) {
                var p = f.p[i*n + j];
                var s = f.m[i*n + j];
                color = getSciColor(p, minP, maxP);
                if (scene.showSmoke) {
                    color[0] = Math.max(0.0, color[0] - 255*s);
                    color[1] = Math.max(0.0, color[1] - 255*s);
                    color[2] = Math.max(0.0, color[2] - 255*s);
                }
            }
            else if (scene.showSmoke) {
                // var s = lightingShader(f.m[i*n+j], [f.u[i*n+j], f.v[i*n+j]], 0.0);
                var s = f.m[i*n + j];
                color[0] = 255*s;
                color[1] = 255*s;
                color[2] = 255*s;
                if (scene.sceneNr == 2)
                    color = getSciColor(s, 0.0, 1.0);
                var s1 = f.m1[i*n + j];
                color1[0] = 255*s1;
                color1[1] = 255*s1;
                color1[2] = 255*s1;
                if(scene.sceneNr == 2){
                    color1 = getSciColor(s1, 0.0, 1.0);
                }
            }
            else if (f.s[i*n + j] == 0.0) {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
                color1[0] = 0;
                color1[1] = 0;
                color1[2] = 0;s
            }

            var x = Math.floor(cX(i * h));
            var y = Math.floor(cY((j+1) * h));
            var cx = Math.floor(cScale * cellScale * h) + 1;
            var cy = Math.floor(cScale * cellScale * h) + 1;

            r = color[0];
            g = color[1];
            b = color[2];
            r1 = color1[0];
            g1 = color1[1];
            b1 = color1[2];

            for (var yi = y; yi < y + cy; yi++) {
                var p = 4 * (yi * canvas.width + x)
                //I think 180 is lowest threshold for flare	and 250 is the lowest threshold for smoke that recreates the smokes in boom beach.
                for (var xi = 0; xi < cx; xi++) {
                    if(r1 < 180  && r < 240){
                        id.data[p++] = (r1 + r)/2;
                        id.data[p++] = (g1 + g)/2;
                        id.data[p++] = (b1 + b)/2;
                        id.data[p++] = 255;
                    } else if(r1 < 180) {
                        id.data[p++] = r1;
                        id.data[p++] = g1;
                        id.data[p++] = b1;
                        id.data[p++] = 255;
                    } else if(r < 240){
                        id.data[p++] = r;
                        id.data[p++] = g;
                        id.data[p++] = b;
                        id.data[p++] = 255;
                    }
                }
            }
        }
    }

    if(scene.realistic_shader){
        id = realisticShader(id, f);
    }

    c.putImageData(id, 0, 0);

    if (scene.showVelocities) {

        c.strokeStyle = "#000000";	
        scale = 0.02;	

        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {

                var u = f.u[i*n + j];
                var v = f.v[i*n + j];

                c.beginPath();

                x0 = cX(i * h);
                x1 = cX(i * h + u * scale);
                y = cY((j + 0.5) * h );

                c.moveTo(x0, y);
                c.lineTo(x1, y);
                c.stroke();

                x = cX((i + 0.5) * h);
                y0 = cY(j * h );
                y1 = cY(j * h + v * scale)

                c.beginPath();
                c.moveTo(x, y0);
                c.lineTo(x, y1);
                c.stroke();

            }
        }
    }

    if (scene.showStreamlines) {

        var segLen = f.h * 0.2;
        var numSegs = 15;

        c.strokeStyle = "#000000";

        for (var i = 1; i < f.numX - 1; i += 5) {
            for (var j = 1; j < f.numY - 1; j += 5) {

                var x = (i + 0.5) * f.h;
                var y = (j + 0.5) * f.h;

                c.beginPath();
                c.moveTo(cX(x), cY(y));

                for (var n = 0; n < numSegs; n++) {
                    var u = f.sampleField(x, y, U_FIELD);
                    var v = f.sampleField(x, y, V_FIELD);
                    l = Math.sqrt(u*u + v*v);
                    // x += u/l * segLen;
                    // y += v/l * segLen;
                    x += u * 0.01;
                    xy += v * 0.01;
                    if (x > f.numX * f.h)
                        break;

                    c.lineTo(cX(x), cY(y));
                }
                c.stroke();
            }
        }
    }

    if (scene.showObstacle) {

        c.strokeW
        r = scene.obstacleRadius + f.h;
        if (scene.showPressure)
            c.fillStyle = "#000000";
        else
            c.fillStyle = "#DDDDDD";
        c.beginPath();	
        c.arc(
            cX(scene.obstacleX), cY(scene.obstacleY), cScale * r, 0.0, 2.0 * Math.PI); 
        c.closePath();
        c.fill();

        c.lineWidth = 3.0;
        c.strokeStyle = "#000000";
        c.beginPath();	
        c.arc(
            cX(scene.obstacleX), cY(scene.obstacleY), cScale * r, 0.0, 2.0 * Math.PI); 
        c.closePath();
        c.stroke();
        c.lineWidth = 1.0;
    }

    if(scene.showGrenade) {
        c.strokeW
        r = scene.grenadeRadius;
    }

    if (scene.showPressure) {
        var s = "pressure: " + minP.toFixed(0) + " - " + maxP.toFixed(0) + " N/m";
        c.fillStyle ="#000000";
        c.font = "16px Arial";
        c.fillText(s, 10, 35);
    }
}

function setObstacle(x, y, reset) {

    var vx = 0.0;
    var vy = 0.0;

    if (!reset) {
        vx = (x - scene.obstacleX) / scene.dt;
        vy = (y - scene.obstacleY) / scene.dt;
    }

    scene.obstacleX = x;
    scene.obstacleY = y;

    var r = scene.obstacleRadius;
    if (!scene.ball) {
        r = 0
    }
    var f = scene.fluid;
    var n = f.numY;
    var cd = Math.sqrt(2) * f.h;

    for (var i = 1; i < f.numX-2; i++) {
        for (var j = 1; j < f.numY-2; j++) {

            f.s[i*n + j] = 1.0;

            dx = (i + 0.5) * f.h - x;
            dy = (j + 0.5) * f.h - y;

            if (dx * dx + dy * dy < r * r) {
                f.s[i*n + j] = 0.0;
                if (scene.sceneNr == 2) {
                    f.m[i*n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr);
                }
                else {
                    f.m[i*n + j] = 1.0;
                }
                f.u[i*n + j] = vx;
                f.u[(i+1)*n + j] = vx;
                f.v[i*n + j] = vy;
                f.v[i*n + j+1] = vy;
            }
        }
    }
    
    if (scene.ball) {
        scene.showObstacle = true;
    } else {
        scene.showObstacle = false; 
    }
}

function setExplosion(x, y) {
   var r = 0.1;
   var f = scene.fluid;
   var n = f.numY;
   var cd = Math.sqrt(2) * f.h;


   for (var i = 1; i < f.numX-2; i++) {
       for (var j = 1; j < f.numY-2; j++) {
           if (true) {
               f.s[i*n + j] = 1.0;


               dx = (i + 0.5) * f.h - x;
               dy = (j + 0.5) * f.h - y;


               if (dx * dx + dy * dy < r * r) {
                   f.s[i*n + j] = 0.0;
                   if (scene.sceneNr == 2)
                       f.m[i*n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                   else
                   if (f.m[i*n + j] < 0.1) {
                       f.m[i*n + j] = 0.1;
                   }
                  
                   // f.u[i*n + j] += 10; //0.9
                   // f.v[i*n + j] += 10; //0.9 added


                   // f.u[(i+1)*n + j] += 10; //0.9
                   // f.v[(i+1)*n + j] += 10; //0.9 added


                   // f.u[i*n + (j+1)] += 10; //0.9
                   // f.v[i*n + (j+1)] += 10; //0.9 added


                   // f.u[(i+1)*n + (j+1)] += 10; //0.9
                   // f.v[(i+1)*n + (j+1)] += 10; //0.9 added


                   var strength = 5;


                   f.u[(i-1)*n + (j-1)] += strength;
                   f.v[(i-1)*n + (j-1)] += strength;


                   f.u[(i-1)*n + (j)] += strength;
                   f.v[(i-1)*n + (j)] += strength;


                   f.u[(i-1)*n + (j+1)] += strength;
                   f.v[(i-1)*n + (j+1)] += strength;


                   f.u[(i)*n + (j-1)] += strength;
                   f.v[(i)*n + (j-1)] += strength;


                   // f.u[(i)*n + (j)] += 10;
                   // f.v[(i)*n + (j)] += 10;


                   f.u[(i+1)*n + (j+1)] += strength;
                   f.v[(i+1)*n + (j+1)] += strength;


                   f.u[(i+1)*n + (j-1)] += strength;
                   f.v[(i+1)*n + (j-1)] += strength;


                   f.u[(i+1)*n + (j)] += strength;
                   f.v[(i+1)*n + (j)] += strength;


                   f.u[(i+1)*n + (j+1)] += strength;
                   f.v[(i+1)*n + (j+1)] += strength;
               }
           }
       }
    }
  }

function setFlareGrenade(x, y) {
    scene.activeFlares.push([x, y, scene.time]);
}

function clearFlareGrenades() {
    var r = 0.01;
    var f = scene.fluid;
    var n = f.numY;
    var cd = Math.sqrt(2) * f.h;
    duration = 3.0;
    var newactiveFlares = []
    var numFlareGrenades = scene.activeFlares.length;
    for (var i = 0; i < numFlareGrenades; i++) { 
        var startTime = scene.activeFlares[i][2];
        if (scene.time - startTime < duration) {
            newactiveFlares.push(scene.activeFlares[i]);
        }
    }
    scene.activeFlares = newactiveFlares;
}

function displayFlareGrenades() {
    var r = scene.flaregrenadeRadius;
       var f = scene.fluid;
       var n = f.numY;
       var cd = Math.sqrt(2) * f.h;
    clearFlareGrenades();
    var numFlareGrenades = scene.activeFlares.length;
    for (var k = 0; k < numFlareGrenades; k++) {
        var currFlare = scene.activeFlares[k];
        var flare_x = currFlare[0];
        var flare_y = currFlare[1];
        for (var i = 1; i < f.numX-2; i++) {
            for (var j = 1; j < f.numY-2; j++) {
                f.s[i*n + j] = 1.0;
                
                dx = (i + 0.5) * f.h - flare_x;
                dy = (j + 0.5) * f.h - flare_y;

                if (dx * dx + dy * dy < r*r) {
                    f.s[i*n + j] = 0.0;
                    if (scene.sceneNr == 2)
                        f.m1[i*n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                    else
                    if (f.m1[i*n + j] > 0.2) {
                        f.m1[i*n + j] = 0.5;
                    }
                    f.u1[i*n + j] = 1/25;
                    f.u1[(i+1)*n + j] = 3/25;
                    f.v1[i*n + j] = 5/25;
                    f.v1[i*n + j+1] = 7/25;
                }
            }
        }
    }
}

function setSmokeGrenade(x, y, time) {
    var r = scene.smokegrenadeRadius;
    var f = scene.fluid;
    var n = f.numY;
    var cd = Math.sqrt(2) * f.h;

    scene.activeSmokes.push([x, y, scene.time]);
    // console.log("scene.activeSmokes: ", scene.activeSmokes);

    for (var i = 1; i < f.numX-2; i++) {
        for (var j = 1; j < f.numY-2; j++) {
            // console.log("f.m[i*n + j]: ", f.m[i*n +j])
            if (true) {
                f.s[i*n + j] = 1.0;

                dx = (i + 0.5) * f.h - x;
                dy = (j + 0.5) * f.h - y;


                if (dx * dx + dy * dy < r * r) {
                    f.s[i*n + j] = 0.0;
                    if (scene.sceneNr == 2)
                        f.m[i*n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                    else
                    if (f.m[i*n + j] > 0.2) {
                        f.m[i*n + j] = 0.90;
                    }
                    f.u[i*n + j] = 0.0001;
                    f.u[(i+1)*n + j] = 0.0001;
                    f.v[i*n + j] = 0.0001;
                    f.v[i*n + j+1] = 0.0001;
                }
                else if (dx * dx + dy * dy < r*r + (0.00325) && dx * dx + dy * dy > r*r - (0.00325)) {
                    f.m[i*n + j] = 1.0;
                }
            }
        }
    }
}

function clearSmokes() {
    var num_activeSmokes = scene.activeSmokes.length;
    // console.log("scene.activeSmokes: ", scene.activeSmokes)
    var removeIndices = []
    for (var i = 0; i < num_activeSmokes; i++) {
        var smoke_x = scene.activeSmokes[i][0];
        var smoke_y = scene.activeSmokes[i][1];
        var smoke_time = scene.activeSmokes[i][2];
        var curr_time = scene.time;

        if (curr_time - smoke_time >= scene.smokegrenadeDuration) {
            clearSmoke(smoke_x, smoke_y);
            removeIndices.push(i);
        }
    }

    var num_removeIndices = removeIndices.length;
    for (var i = 0; i < num_removeIndices; i++) {
        scene.activeSmokes.splice(removeIndices[i], 1);
    }
}

function clearSmoke(x, y) {
    var r = scene.smokegrenadeRadius;
    var f = scene.fluid;
    var n = f.numY;
    var cd = Math.sqrt(2) * f.h;

    for (var i = 1; i < f.numX-2; i++) {
        for (var j = 1; j < f.numY-2; j++) {
            f.s[i*n + j] = 1.0;

            var clear = true; 

            var num_activeSmokes = scene.activeSmokes.length;
            if(num_activeSmokes == 1){
                f.m[i*n + j] = 1.0;
            }
            for (var k = 0; k < num_activeSmokes; k++) {
                var smoke_x = scene.activeSmokes[k][0];
                var smoke_y = scene.activeSmokes[k][1];

                dx = (i + 0.5) * f.h - smoke_x;
                dy = (j + 0.5) * f.h - smoke_y;

                var store = dx * dx + dy * dy > r*r +(0.004);
                clear = clear && store;
            }

            if (clear) {
                f.m[i*n + j] = 1.0;
            }
        }
    }
}



// interaction -------------------------------------------------------

var mouseDown = false;

function startDrag(x, y) {
    let bounds = canvas.getBoundingClientRect();

    let mx = x - bounds.left - canvas.clientLeft;
    let my = y - bounds.top - canvas.clientTop;
    mouseDown = true;

    x = mx / cScale;
    y = (canvas.height - my) / cScale;

    setObstacle(x,y, true);
}

function startExplosion(x, y) { //added
    let bounds = canvas.getBoundingClientRect();

    let mx = x - bounds.left - canvas.clientLeft;
    let my = y - bounds.top - canvas.clientTop;
    mouseDown = true; 

    x = mx / cScale;
    y = (canvas.height - my) / cScale; 

    setExplosion(x, y);
}

function startFlareGrenade(x, y) {
    let bounds = canvas.getBoundingClientRect();

    let mx = x - bounds.left - canvas.clientLeft;
    let my = y - bounds.top - canvas.clientTop;
    mouseDown = true; 

    x = mx / cScale;
    y = (canvas.height - my) / cScale; 

    setFlareGrenade(x, y);
}

function startSmokeGrenade(x, y) {
    let bounds = canvas.getBoundingClientRect();

    let mx = x - bounds.left - canvas.clientLeft;
    let my = y - bounds.top - canvas.clientTop;
    mouseDown = true; 

    x = mx / cScale;
    y = (canvas.height - my) / cScale; 

    setSmokeGrenade(x, y);
}

function drag(x, y) {
    if (mouseDown) {
        let bounds = canvas.getBoundingClientRect();
        let mx = x - bounds.left - canvas.clientLeft;
        let my = y - bounds.top - canvas.clientTop;
        x = mx / cScale;
        y = (canvas.height - my) / cScale;
        if (scene.ball) {
            setObstacle(x,y, false);
        }
    }
}

function endDrag() {
    mouseDown = false;
}

canvas.addEventListener('mousedown', event => {
    if (scene.ball) {
        startDrag(event.x, event.y);
    }
    else if (scene.explosion) {
        startExplosion(event.x, event.y);
    }
    else if (scene.flaregrenade){
        startFlareGrenade(event.x, event.y);
    }
    else if (scene.smokegrenade) {
        startSmokeGrenade(event.x, event.y);
    }
});

canvas.addEventListener('mouseup', event => {
    endDrag();
});

canvas.addEventListener('mousemove', event => {
    drag(event.x, event.y);
});

canvas.addEventListener('touchstart', event => {
    startDrag(event.touches[0].clientX, event.touches[0].clientY)
});

canvas.addEventListener('touchend', event => {
    endDrag()
});

canvas.addEventListener('touchmove', event => {
    event.preventDefault();
    event.stopImmediatePropagation();
    drag(event.touches[0].clientX, event.touches[0].clientY)
}, { passive: false});


document.addEventListener('keydown', event => {
    switch(event.key) {
        case 'p': scene.paused = !scene.paused; break;
        case 'm': scene.paused = false; simulate(); scene.paused = true; break;
    }
});

function toggleStart()
{
    var button = document.getElementById('startButton');
    if (scene.paused)
        button.innerHTML = "Stop";
    else
        button.innerHTML = "Start";
    scene.paused = !scene.paused;
}

// main -------------------------------------------------------

function simulate() 
{
    if (!scene.paused) {
        scene.time += scene.dt;
        scene.fluid.simulate(scene.dt, scene.gravity, scene.numIters);
        scene.frameNr++;
    }

    // if(scene.activeFlares.length > 0){
    // 	for(var i = 0; i < scene.activeFlares.length; i++){
    // 		setFlareGrenade(scene.activeFlares[i][0], scene.activeFlares[i][1]);
    // 		var curr_time = scene.time;
    // 		var duration = 2.0;
    // 		if (curr_time - scene.activeFlares[i][2] >= duration) {
    // 			scene.activeFlares.splice(i, 1);
    // 		}
    // 	}
    // }

    if (scene.activeFlares.length > 0) {
        displayFlareGrenades();
    }

    if (!scene.ball) {
        setObstacle(0.4, 0.5, true)
    }
    if(scene.activeSmokes.length > 0){
        clearSmokes();
    }
}

function update() {
    simulate();
    draw();
    requestAnimationFrame(update);
}

setupScene(4);
    update();