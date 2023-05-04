/* When the user clicks on the button,
toggle between hiding and showing the dropdown content */
function myFunction() {
    document.getElementById("myDropdown").classList.toggle("show");
}

function myFunction2() {
    document.getElementById("myDropdown2").classList.toggle("show");
}

// Close the dropdown if the user clicks outside of it
window.onclick = function (event) {
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

var simHeight = 1.1;
var cScale = canvas.height / simHeight;
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

class Fluid {
    constructor(density, numX, numY, h) {
        this.baseDensity = density;
        this.density = new Float32Array(this.numCells);
        this.density.fill(density);

        this.numX = numX + 2;
        this.numY = numY + 2;
        this.numCells = this.numX * this.numY;
        this.h = h;
        this.u = new Float32Array(this.numCells);
        this.v = new Float32Array(this.numCells);
        this.w = Array(this.numCells).fill([0, 0]);
        this.newU = new Float32Array(this.numCells);
        this.newV = new Float32Array(this.numCells);
        this.newW = new Array(this.numCells).fill([0, 0]);
        this.normal = Array(this.numCells).fill([0, 0]);
        this.force = Array(this.numCells).fill([0, 0]);
        this.p = new Float32Array(this.numCells);
        this.s = new Float32Array(this.numCells);
        this.m = new Float32Array(this.numCells);
        this.newM = new Float32Array(this.numCells);
        this.m.fill(1.0)
        var num = numX * numY;

        this.temp = new Float32Array(this.numCells);
        this.temp.fill(0.3); // 0.3 is smoke, 0.6 is fire, 1 would be an orange hot fire (this is vaguely mimicking the real temps in c)
    }

    updateDensity() {
        // density is initially set to 1000 (I think this is an arbitrary number they chose)
        // density = p / (R * T)
        var n = this.numY;
        for (var i = 1; i < this.numX - 1; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
                // using 0.3 as a baseline, any difference from that will cause a difference in density
                // higher temp means lower density (this is technically what causes buoyancy)
                this.density[i * n + j] += ((- (this.temp[i * n + j] - 0.3) / 0.3) * 0.5) * this.baseDensity;
            }
        }
    }

    updateTemp() {
        // T = p / (R * density)
        var n = this.numY;
        for (var i = 1; i <= this.numX - 1; i++) {
            for (var j = 1; j <= this.numY - 1; j++) {
                if (this.p[i * n + j] != 0) {
                    // radiation
                    var convConst = 0.1;
                    if (this.u[i * n + j] > 0) {
                        this.temp[(i * n) + j] += convConst * this.u[i * n + j] * (this.temp[((i - 1) * n) + j] - this.temp[(i * n) + j])
                    } else {
                        this.temp[(i * n) + j] += convConst * this.u[i * n + j] * (this.temp[((i + 1) * n) + j] - this.temp[(i * n) + j])
                    }
                    if (this.v[i * n + j] > 0) {
                        this.temp[(i * n) + j] += convConst * this.v[i * n + j] * (this.temp[(i * n) + j - 1] - this.temp[(i * n) + j])
                    } else {
                        this.temp[(i * n) + j] += convConst * this.v[i * n + j] * (this.temp[(i * n) + j + 1] - this.temp[(i * n) + j])
                    }

                    var upMovement = 0.3;
                    var downMovement = 0.1;
                    var leftMovement = 0.2;
                    var rightMovement = 0.16;
                    this.temp[(i * n) + j] += downMovement * (this.temp[(i * n) + j + 1] - this.temp[(i * n) + j])
                    this.temp[(i * n) + j] += upMovement * (this.temp[(i * n) + j - 1] - this.temp[(i * n) + j])
                    this.temp[(i * n) + j] += leftMovement * (this.temp[((i + 1) * n) + j] - this.temp[(i * n) + j])
                    this.temp[(i * n) + j] += rightMovement * (this.temp[((i - 1) * n) + j] - this.temp[(i * n) + j])

                    if (this.temp[(i * n) + j] > 1.0) {
                        this.temp[(i * n) + j] = 1.0;
                    }
                }
            }
        }
    }

    addHeat() {
        var n = this.numY;
        for (var i = 1; i < this.numX - 2; i++) {
            for (var j = 1; j < this.numY - 2; j++) {
                if (this.s[i * n + j] == 0) {
                    for (var l = -1; l <= 1; l++) {
                        for (var m = -1; m <= 1; m++) {
                            if ((i + l) < this.numX - 2 && ((i + l) > 1)) {
                                if (j + m < this.numY - 2 && (j + m > 1)) {
                                    this.temp[(i + l) * n + j + m] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    addSmoke() {
        var n = this.numY;
        for (var i = 1; i < this.numX - 2; i++) {
            for (var j = 1; j < this.numY - 2; j++) {
                if (this.s[i * n + j] == 0) {
                    for (var l = -1; l <= 1; l++) {
                        for (var m = -1; m <= 1; m++) {
                            if ((i + l) < this.numX - 2 && ((i + l) > 1)) {
                                if (j + m < this.numY - 2 && (j + m > 1)) {
                                    if (this.s[(i + l) * n + j + m] == 1) {
                                        this.m[(i + l) * n + j + m] = 0.2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    tempReductionOverTime() {
        var n = this.numY;
        var reductionConst = 0.02;
        for (var i = 1; i < this.numX - 1; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
                this.temp[i * n + j] += (0.3 - this.temp[i * n + j]) * reductionConst;
            }
        }
    }

    integrate(dt, gravity) {
        var n = this.numY;
        var buoyancyConst = 0.75;
        for (var i = 1; i < this.numX; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
                if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0) {
                    this.v[i * n + j] += gravity * dt;
                    if (this.temp[(i * n) + j] - this.temp[(i * n) + j + 1] > 0) {
                        this.v[i * n + j] += (buoyancyConst * (this.temp[(i * n) + j] - this.temp[(i * n) + j + 1]));
                    }
                }
            }
        }
    }

    solveIncompressibility(numIters, dt) {

        var n = this.numY;
        for (var iter = 0; iter < numIters; iter++) {
            for (var i = 1; i < this.numX - 1; i++) {
                for (var j = 1; j < this.numY - 1; j++) {
                    // if (this.s[i*n + j] == 0.0)
                    // 	continue;
                    if (this.s[i * n + j] != 0.0) {
                        var s = this.s[i * n + j];
                        var sx0 = this.s[(i - 1) * n + j];
                        var sx1 = this.s[(i + 1) * n + j];
                        var sy0 = this.s[i * n + j - 1];
                        var sy1 = this.s[i * n + j + 1];
                        var s = sx0 + sx1 + sy0 + sy1;
                        // if (s == 0.0)
                        // 	continue;
                        if (s != 0.0) {
                            var div = this.u[(i + 1) * n + j] - this.u[i * n + j] +
                                this.v[i * n + j + 1] - this.v[i * n + j];

                            var cp = this.density[i * n + j] * this.h / dt;
                            var p = -div / s;
                            p *= scene.overRelaxation;

                            this.p[i * n + j] += cp * p;

                            this.u[i * n + j] -= sx0 * p;
                            this.u[(i + 1) * n + j] += sx1 * p;
                            this.v[i * n + j] -= sy0 * p;
                            this.v[i * n + j + 1] += sy1 * p;
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
        for (var i = 1; i < this.numX - 1; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
                // if (this.s[i*n + j] == 0.0)
                // 	continue;
                if (this.grenadeX == i && this.grenadeY == j) {
                    this.p[i * n + j] = 10;
                }
            }
        }
    }

    extrapolate() {
        var n = this.numY;
        for (var i = 0; i < this.numX; i++) {
            this.u[i * n + 0] = this.u[i * n + 1];
            this.u[i * n + this.numY - 1] = this.u[i * n + this.numY - 2];
        }
        for (var j = 0; j < this.numY; j++) {
            this.v[0 * n + j] = this.v[1 * n + j];
            this.v[(this.numX - 1) * n + j] = this.v[(this.numX - 2) * n + j]
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

        switch (field) {
            case U_FIELD: f = this.u; dy = h2; break;
            case V_FIELD: f = this.v; dx = h2; break;
            case S_FIELD: f = this.m; dx = h2; dy = h2; break

        }

        var x0 = Math.min(Math.floor((x - dx) * h1), this.numX - 1);
        var tx = ((x - dx) - x0 * h) * h1;
        var x1 = Math.min(x0 + 1, this.numX - 1);

        var y0 = Math.min(Math.floor((y - dy) * h1), this.numY - 1);
        var ty = ((y - dy) - y0 * h) * h1;
        var y1 = Math.min(y0 + 1, this.numY - 1);

        var sx = 1.0 - tx;
        var sy = 1.0 - ty;

        var val = sx * sy * f[x0 * n + y0] +
            tx * sy * f[x1 * n + y0] +
            tx * ty * f[x1 * n + y1] +
            sx * ty * f[x0 * n + y1];

        return val;
    }

    avgU(i, j) {
        var n = this.numY;
        var u = (this.u[i * n + j - 1] + this.u[i * n + j] +
            this.u[(i + 1) * n + j - 1] + this.u[(i + 1) * n + j]) * 0.25;
        return u;
    }

    avgV(i, j) {
        var n = this.numY;
        var v = (this.v[(i - 1) * n + j] + this.v[i * n + j] +
            this.v[(i - 1) * n + j + 1] + this.v[i * n + j + 1]) * 0.25;
        return v;
    }

    advectVel(dt) {

        this.newU.set(this.u);
        this.newV.set(this.v);

        var n = this.numY;
        var h = this.h;
        var h2 = 0.5 * h;

        for (var i = 1; i < this.numX; i++) {
            for (var j = 1; j < this.numY; j++) {

                cnt++;

                // u component
                if (this.s[i * n + j] != 0.0 && this.s[(i - 1) * n + j] != 0.0 && j < this.numY - 1) {
                    var x = i * h;
                    var y = j * h + h2;
                    var u = this.u[i * n + j];
                    var v = this.avgV(i, j);
                    //						var v = this.sampleField(x,y, V_FIELD);
                    x = x - dt * u;
                    y = y - dt * v;
                    u = this.sampleField(x, y, U_FIELD);
                    this.newU[i * n + j] = u;
                }
                // v component
                if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0 && i < this.numX - 1) {
                    var x = i * h + h2;
                    var y = j * h;
                    var u = this.avgU(i, j);
                    //						var u = this.sampleField(x,y, U_FIELD);
                    var v = this.v[i * n + j];
                    x = x - dt * u;
                    y = y - dt * v;
                    v = this.sampleField(x, y, V_FIELD);
                    this.newV[i * n + j] = v;
                }
            }
        }

        this.u.set(this.newU);
        this.v.set(this.newV);
    }

    advectSmoke(dt) {

        this.newM.set(this.m);

        var n = this.numY;
        var h = this.h;
        var h2 = 0.5 * h;

        for (var i = 1; i < this.numX - 1; i++) {
            for (var j = 1; j < this.numY - 1; j++) {

                if (this.s[i * n + j] != 0.0) {
                    var u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5;
                    var v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5;
                    var x = i * h + h2 - dt * u;
                    var y = j * h + h2 - dt * v;

                    this.newM[i * n + j] = this.sampleField(x, y, S_FIELD);
                }
            }
        }
        this.m.set(this.newM);
    }

    curl(x, y, direction = 1) {
        var n = this.numY;
        var my_curl = [0, 0];
        my_curl[0] = -this.u[x * n + y + 1] + this.u[x * n + y - 1]; //counterclockwise rotation
        my_curl[1] = -this.v[(x - 1) * n + y] + this.v[(x + 1) * n + y]; //counterclockwise rotation
        // console.log(my_curl);
        if (direction != 1) { //clockwise rotation
            my_curl[0] = -my_curl[0];
            my_curl[1] = -my_curl[1];
        }
        return my_curl;
    }

    norm(arr) { //l-2 norm of array/vector
        var sum = 0;
        for (var i = 0; i < arr.length; i++) {
            sum += Math.pow(arr[i], 2);
        }
        return Math.sqrt(sum);
    }

    dot(arr1, arr2) { //dot product of two vectors
        var out = 0;
        for (var i = 0; i < arr1.length; i++) {
            out += arr1[i] * arr2[i];
        }
        return out;
    }

    cross(arr1, arr2) { //magnitude of cross product of two vectors
        // console.log(this.dot(arr1, arr2) / (this.norm(arr1) * this.norm(arr2)));
        return this.norm(arr1) * this.norm(arr2) * Math.sin(Math.acos(this.dot(arr1, arr2) / (this.norm(arr1) * this.norm(arr2))));
    }

    vorticityConfinement(dt) {
        
        const cartesian = (arr1, arr2) => {
        const res = [];
        for(let i = 0; i < arr1.length; i++){
            for(let j = 0; j < arr2.length; j++){
                res.push(
                    [arr1[i]].concat(arr2[j])
                );
            };
        };
        return res;
        };

        this.newU.set(this.u);
        this.newV.set(this.v);
        for (var i = 0; i < this.w.length; i++) {
            this.newW[i] = this.w[i];
        }
        // this.newW.set(this.w);

        // var f = scene.fluid;

        var width = this.numX;
        var height = this.numY;

        var dx = 0.0;
        var dy = 0.0;
        var len = 0.0;
        var vorticity = 1.3;
        
        var abs_grad = [0, 0];
        var my_norm = 0;

        var n = this.numY;

        for (var my_y = 2; my_y < this.numY - 2; my_y++) {
            for (var my_x = 2; my_x < this.numX - 2; my_x++) {
                this.newW[my_x * n + my_y] = this.curl(my_x, my_y);
                for (var i = 0; i < abs_grad.length; i++) {
                    // console.log(this.newW[my_x * n + my_y][i] - this.w[my_x * n + my_y][i]);
                    abs_grad[i] = Math.abs(this.newW[my_x * n + my_y][i] - this.w[my_x * n + my_y][i]) / dt;
                    // console.log(abs_grad[i]);
                }
                my_norm = this.norm(abs_grad);
                for (var i = 0; i < this.normal[my_x * n + my_y].length; i++) {
                    // console.log(abs_grad[i]);
                    this.normal[my_x * n + my_y][i] = abs_grad[i] / my_norm;
                }
                if (this.normal[my_x * n + my_y][0] == 0 && this.normal[my_x * n + my_y][1] == 0) {
                    break;
                }
                // console.log(this.normal[my_x * n + my_y]);
                for (var i = 0; i < this.force[my_x * n + my_y].length; i++) {
                    this.force[my_x * n + my_y][i] = vorticity * this.h * this.cross(this.normal[my_x * n + my_y], this.newW[my_x * n + my_y]);
                }
                // console.log(this.cross(this.normal[my_x * n + my_y], this.newW[my_x * n + my_y]));
            }
        }

        var curr_mass = 0;
        var s = 0;
        var sx0 = 0;
        var sx1 = 0;
        var sy0 = 0;
        var sy1 = 0;
        var p = 0;
        var epsilon = Math.pow(10, -4);
        var gamma = Math.pow(10, 0);
        for (var my_y = 2; my_y < this.numY - 2; my_y++) {
            for (var my_x = 2; my_x < this.numX - 2; my_x++) {
                if (this.s[my_x * n + my_y] != 0.0) {
                    s = this.s[my_x * n + my_y];
                    sx0 = this.s[(my_x - 1) * n + my_y];
                    sx1 = this.s[(my_x + 1) * n + my_y];
                    sy0 = this.s[my_x * n + (my_y - 1)];
                    sy1 = this.s[my_x * n + (my_y + 1)];
                    s = sx0 + sx1 + sy0 + sy1;
                    // if (s == 0.0)
                    // 	continue;
                    if (s != 0.0) {
                        var div = this.u[(my_x + 1) * n + my_y] - this.u[my_x * n + my_y] + 
                            this.v[my_x * n + (my_y + 1)] - this.v[my_x * n + my_y];

                        p = -div / s;
                        p *= scene.overRelaxation;
                        // this.p[my_x * n + my_y] += cp * p;

                        // this.u[my_x * n + my_y] -= sx0 * p;
                        // this.u[(my_x + 1) * n + my_y] += sx1 * p;
                        // this.v[my_x * n + my_y] -= sy0 * p;
                        // this.v[my_x * n + my_y + 1] += sy1 * p;
                    }
                }
                curr_mass = 0;
                if (sx0 * p > 0) {
                    curr_mass += sx0 * p;
                } else {
                    curr_mass += sx1 * p;
                }
                if (sy0 * p > 0) {
                    curr_mass += sy0 * p;
                } else {
                    curr_mass += sy1 * p;
                }
                // curr_mass = (sx0 * p + sy0 * p);
                // if (curr_mass != 0) {
                // 	console.log(curr_mass);
                // }
                if (!isNaN(curr_mass) && curr_mass > epsilon && curr_mass < gamma && this.force[my_x * n + my_y][0] < gamma && this.force[my_x * n + my_y][1] < gamma) {
                    // console.log("force" + this.force[my_x * n + my_y][0].toString());
                    // console.log(f.m[my_x * n + my_y]);
                    if (this.u[my_x * n + my_y] > 0 && this.force[my_x * n + my_y][0] > 0 || this.u[my_x * n + my_y] < 0 && this.force[my_x * n + my_y][0] < 0) {
                        this.newU[my_x * n + my_y] = this.u[my_x * n + my_y] + vorticity * dt * this.force[my_x * n + my_y][0] / curr_mass;
                    }
                    // console.log(f.m[my_x * n + my_y]);
                    if (this.v[my_x * n + my_y] > 0 && this.force[my_x * n + my_y][1] > 0 || this.v[my_x * n + my_y] < 0 && this.force[my_x * n + my_y][1] < 0) {
                        this.newV[my_x * n + my_y] = this.v[my_x * n + my_y] + vorticity * dt * this.force[my_x * n + my_y][1] / curr_mass;
                    }
                }
            }
        }

        this.u.set(this.newU);
        this.v.set(this.newV);
        for (var i = 0; i < this.w.length; i++) {
            this.w[i] = this.newW[i];
        }
    }

    // ----------------- end of simulator ------------------------------


    simulate(dt, gravity, numIters) {
        this.tempReductionOverTime();
        this.addHeat();
        this.addSmoke();

        if (this.showGrenade) {
            this.integrate(dt, gravity);
            this.p.fill(0, 0);
            this.grenadeCompressibility(x, y);
            this.extrapolate();
            this.advectVel(dt);
            this.vorticityConfinement(dt);
            this.advectSmoke(dt);
        }

        this.integrate(dt, gravity);

        this.updateTemp();
        this.updateDensity();

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
    ggravity: -9.81,
    dt: 1.0 / 120.0,
    numIters: 100,
    frameNr: 0,
    overRelaxation: 1.9,
    obstacleX: 0.0,
    obstacleY: 0.0,
    obstacleRadius: 0.15,
    paused: false,
    sceneNr: 0,
    showObstacle: false,
    showStreamlines: false,
    showVelocities: false,
    showPressure: false,
    showTemp: false,
    showSmoke: true,
    showGrenade: false,
    fire: false,
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

function setupScene(sceneNr = 0) {
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

    f = scene.fluid = new Fluid(density, numX, numY, h);

    var n = f.numY;

    if (sceneNr == 0) {   		// tank

        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {
                var s = 1.0;	// fluid
                if (i == 0 || i == f.numX - 1 || j == 0)
                    s = 0.0;	// solid
                f.s[i * n + j] = s
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
        scene.blocky = false;
        scene.realistic = false;
        scene.color = false;
    }
    else if (sceneNr == 1 || sceneNr == 3) { // vortex shedding
        var inVel = 2.0;
        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {
                var s = 1.0;	// fluid
                if (i == 0 || j == 0 || j == f.numY - 1)
                    s = 0.0;	// solid
                f.s[i * n + j] = s

                if (i == 1) {
                    f.u[i * n + j] = inVel;
                }
            }
        }

        var pipeH = 0.1 * f.numY;
        var minJ = Math.floor(0.5 * f.numY - 0.5 * pipeH);
        var maxJ = Math.floor(0.5 * f.numY + 0.5 * pipeH);

        for (var j = minJ; j < maxJ; j++)
            f.m[j] = 0.0;


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
        scene.blocky = false;
        scene.realistic = false;
        scene.color = false;
        scene.explosion = false; //added
        scene.flaregrenade = false; //added
        scene.smokegrenade = false; //added

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
        scene.blocky = false;
        scene.realistic = false;
        scene.color = false;
        scene.explosion = false; //added
        scene.flaregrenade = false; //added
        scene.smokegrenade = false; //added

    } else if (sceneNr == 4) { // None
        scene.gravity = 0.0;
        scene.overRelaxation = 1.0;
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;
        scene.obstacleRadius = 0.1;
        scene.ball = false;
        scene.showTemp = false;
        scene.activeFlares = []
        scene.explosion = false; //added
        scene.flaregrenade = false; //added
        scene.smokegrenade = false; //added

        scene.smokeAppOne = true; // added
        scene.smokeAppTwo = false; // added
        scene.smokeAppThree = false; // added
        scene.fire = false;
        scene.blocky = false;
        scene.realistic = false;
        scene.color = false;
    }

    document.getElementById("smokeButton").checked = scene.showSmoke;
    document.getElementById("overrelaxButton").checked = scene.overRelaxation > 1.0;
    document.getElementById("ballButton").checked = scene.ball;
    document.getElementById("flaregrenadeButton").checked = scene.flaregrenade;
    document.getElementById("smokegrenadeButton").checked = scene.smokegrenade;
    document.getElementById("tempButton").checked = scene.showTemp;
    document.getElementById("smokeAppOne").checked = scene.smokeAppOne;
    document.getElementById("smokeAppTwo").checked = scene.smokeAppTwo;
    document.getElementById("smokeAppThree").checked = scene.smokeAppThree;
    document.getElementById("fire").checked = scene.fire;
    document.getElementById("blocky").checked = scene.blocky;
    document.getElementById("color").checked = scene.color;
}


// draw -------------------------------------------------------

function setColor(r, g, b) {
    c.fillStyle = `rgb(
        ${Math.floor(255 * r)},
        ${Math.floor(255 * g)},
        ${Math.floor(255 * b)})`
    c.strokeStyle = `rgb(
        ${Math.floor(255 * r)},
        ${Math.floor(255 * g)},
        ${Math.floor(255 * b)})`
}

function getSciColor(val, minVal, maxVal) {
    val = Math.min(Math.max(val, minVal), maxVal - 0.0001);
    var d = maxVal - minVal;
    val = d == 0.0 ? 0.5 : (val - minVal) / d;
    var m = 0.25;
    var num = Math.floor(val / m);
    var s = (val - num * m) / m;
    var r, g, b;

    switch (num) {
        case 0: r = 0.0; g = s; b = 1.0; break;
        case 1: r = 0.0; g = 1.0; b = 1.0 - s; break;
        case 2: r = s; g = 1.0; b = 0.0; break;
        case 3: r = 1.0; g = 1.0 - s; b = 0.0; break;
    }

    return [255 * r, 255 * g, 255 * b, 255]
}

function norm(vector) {
    var vectorLength = vector.length;
    var total = 0.0;
    for (var i = 0; i < vectorLength; i++) {
        total += vector[i] * vector[i];
    }
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

function vecMul(vector1, vector2) {
    final = []
    var vectorLength = vector1.length;
    for (var i = 0; i < vectorLength; i++) {
        final.push(vector1[i] * vector2[i]);
    }
    return final
}

function vecscalarMul(vector, scalar) {
    final = []
    var vectorLength = vector.length;
    for (var i = 0; i < vectorLength; i++) {
        final.push(vector[i] * scalar);
    }
    return final
}

function vecAdd(vectorofVectors) {
    final = []
    var numVectors = vectorofVectors.length;
    var numElements = vectorofVectors[0].length;

    for (var i = 0; i < numElements; i++) {
        store = 0;
        for (var j = 0; j < numVectors; j++) {
            store += vectorofVectors[j][i];
        }
        final.push(store);
    }

    return final;
}

function vecSub(vector1, vector2) {
    final = [];
    var numElements = vector1.length;

    for (var i = 0; i < numElements; i++) {
        final.push(vector1[i] - vector2[i]);
    }
    return final;
}

// function lightingShader(mass, velocity, ambient = 0.2) {
//     var normVelocity = normalize(velocity);
//     var normLightDir = normalize([1.0, 1.0]);
//     var store = dot(normVelocity, normLightDir);
//     // console.log("store: ", store);

//     var shading = Math.max(store, 0.0) + ambient;
//     shading = Math.min(shading, 1.0);
//     shading = mass;
//     var colorValue = [255 * shading, 255 * shading, 255 * shading];
//     return colorValue;
// }

function vangoghShader(mass) { 
    if (mass > 0.95) {
        //Black
        color = [0, 0, 0, mass * 255];
    } 
    else if (mass >= 0.6) {
        //Blueish Portion
        if (mass >= 0.6 && mass < 0.7) {
            //Light Blue
            color = [0, 0, 255, 150]; //150
        }
        else if (mass >= 0.7 && mass < 0.77) {
            //Black
            color = [0, 0, 0, 255]; //255 
        }
        else if (mass >= 0.77 && mass < 0.85) {
            //White
            color = [255, 255, 255, 255]; //255
        }
        else if (mass >= 0.85 && mass < 0.90) {
            //Blueish Purple
            color = [115, 155, 255, 255]; //255
        }
        else {
            // //Dark Blue
            color = [0, 150, 255, 255]; //255
            // // Blueish Purple
            // color = [115, 155, 255, 255]; //255
        }
    }
    else if (mass >= 0.3) {
        //Yellowish
        if (mass >= 0.5 && mass < 0.6) {
            // //Greenish yellow
            // color = [200, 255, 0, 127]
            //Blueish Purple
            color = [115, 155, 255, 255]; //255
        }
        else if (mass >= 0.4 && mass < 0.5) {
            // //Lighter Yellow
            // color = [255, 255, 0, 127]
            //Blueish Purple but 1
            color = [115, 155, 255, 255]; //255
        }
        else {
            //Light Yellow
            color = [255, 255, 0, 255]; //255
        }
        
    }
    else {
        //Cyanish
        if (mass >= 0.2 && mass < 0.3) {
            // //Lightish blue
            // color = [0, 155, 200, 200];
            // //White
            // color = [255, 255, 255, 255]
            //Blueish Purple
            color = [115, 155, 255, 255]; //255
        }
        else if (mass >= 0.1 && mass < 0.2) {
            // //Purpleish
            // color = [50, 155, 200, 200]
            //Blueish Purple
            color = [115, 155, 255, 255]; //255
        }
        else {
            //Cyan
            color = [0, 190, 255, 255]; //255
        }
    }
    return color;
}

function lerp(t, color1, color2) {
    newColor = [(t * color1[0]) + (1 - t) * color2[0], (t * color1[1]) + (1 - t) * color2[1], (t * color1[2]) + (1 - t) * color2[2], 255];
    return newColor;
}

function fireShader(heat) {
    // Dark Orange
    color1 = [152, 34, 7, 255];

    // lighter orange
    color2 = [222, 64, 6, 255];

    // orange
    color3 = [245, 100, 8, 255];

    // yellow 
    color4 = [249, 187, 11, 255];

    // light orange
    color5 = [248, 137, 9, 255];

    // very pale yellow
    color6 = [250, 250, 163, 255];

    if (heat >= 0.95) {
        // Black 
        color = [0, 0, 0, 255];
    }
    else if (heat >= 0.9) {
        color = lerp((heat - 0.9) / 0.05, color, color1);
    }
    else if (heat >= 0.85) {
        color = lerp((heat - 0.85) / 0.05, color1, color2);
    }
    else if (heat >= 0.75) {
        color = lerp((heat - 0.75) / 0.1, color2, color3);
    }
    else if (heat >= 0.65) {
        color = lerp((heat - 0.65) / 0.1, color3, color4);
    }
    else if (heat >= 0.55) {
        color = lerp((heat - 0.55) / 0.1, color4, color5);
    }
    else if (heat >= 0.45) {
        color = lerp((heat - 0.45) / 0.1, color5, color6);
    }
    else {
        color = [255, 255, 255, 0];
    }
    return color;
}

function blockyShader(mass, smokeColor) {
    if (mass > 0.9) {
        color = [255, 255, 255, 0];
    } else {
        color = [0, 0, 0, 255];
    }
    return color;
}

function realisticShader(id, f) {
    var cx = Math.floor(cScale * 1.1 * f.h) + 1;
    var cy = Math.floor(cScale * 1.1 * f.h) + 1;
    var copy = id;
    for (var i = 0; i < f.numX; i++) {
        for (var j = 0; j < f.numY; j++) {
            var x = Math.floor(cX(i * f.h));
            var y = Math.floor(cY((j + 1) * f.h));
            for (var yi = y; yi < y + cy; yi++) {
                var p = 4 * (yi * canvas.width + x);
                for (var xi = 0; xi < cx; xi++) {
                    if (p - 3 >= 0 && id.length > p + 4 * canvas.width + 7 && p - 4 * canvas.width - 3 >= 0) {
                        avg_r = id[p + 1] + id[p - 3] + id[p + 5] + id[p + 4(canvas.width) + 1] + id[p - 4 * canvas.width + 1]+ id[p + 4 * (canvas.width) + 5] + id[p + 4 * (canvas.width) - 3] + id[p - 4 * canvas.width + 5]+ + id[p - 4 * canvas.width - 3];
                        avg_g = id[p + 2] + id[p - 2] + id[p + 6] + id[p + 4(canvas.width) + 2] + id[p - 4 * canvas.width + 2]+ id[p + 4 * (canvas.width) + 6] + id[p + 4 * (canvas.width) - 2] + id[p - 4 * canvas.width + 6]+ + id[p - 4 * canvas.width - 2];
                        avg_b = id[p + 3] + id[p - 1] + id[p + 7] + id[p + 4(canvas.width) + 3] + id[p - 4 * canvas.width + 3]+ id[p + 4 * (canvas.width) + 7] + id[p + 4 * (canvas.width) - 1] + id[p - 4 * canvas.width + 7]+ + id[p - 4 * canvas.width - 1];
                        avg_a =
                            copy[p++] = avg_r / 9;
                        copy[p++] = avg_g / 9;
                        copy[p++] = avg_b / 9;
                        copy[p++] = id[p + 4];
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

function colorShader(mass, smokeColor) {
    // if (mass > 0.97) {
    //     color = [255, 255, 255, 255];
    // }
    // else {
    //     color = [(1.0 - mass) * smokeColor[0], (1. - mass)smokeColor[1], (1.0 - mass)smokeColor[2], (1.0 - mass)smokeColor[3]]
    // }
    color = [(1.0 - mass) * smokeColor[0], (1.0 - mass) * smokeColor[1], (1.0 - mass) * smokeColor[2], 255];
    return color;
}

function colorShader(mass, smokeColor) {
    // if (mass > 0.97) {
    //     color = [255, 255, 255, 255];
    // }
    // else {
    //     color = [(1.0 - mass) * smokeColor[0], (1. - mass)*smokeColor[1], (1.0 - mass)*smokeColor[2], (1.0 - mass)*smokeColor[3]]
    // }
    color = [(1.0 - mass) * smokeColor[0], (1.0 - mass)*smokeColor[1], (1.0 - mass)*smokeColor[2], (1.0 - mass)*smokeColor[3]]
    return color;
}

function phongShader(mass, velocity, loc) {
    // console.log("mass: ", mass);
    // console.log("velocity: ", velocity);
    // console.log("log: ", loc);
    var lightPos = [0.0, 0.0, 0.0];
    var camPos = [scene.fluid.numX / 2.0, scene.fluid.numY / 2.0, 10.0];

    var ka = 0.1;
    var kd = [1.0, 1.0, 1.0];
    var ks = 0.5;
    var Ia = [1.0, 1.0, 1.0];
    var p = 100;

    //Ambient Lighting
    var ambLighting = vecscalarMul(Ia, ka);

    //Diffuse Reflection
    var v_normal = [loc[0], loc[1], 1.0]
    var v_normal_normalize = normalize(v_normal); //The v normal is just pointing up
    var lightDis = norm(vecSub(loc, lightPos));
    var lightVec = normalize(vecSub(loc, lightPos));
    // console.log("v_normal_normalize: ", v_normal_normalize);
    // console.log("lightVec: ", lightVec)
    var m = Math.max(0.0, dot(v_normal_normalize, lightVec));
    // console.log("m: ", m);
    if (lightDis == 0) {
        var I_r = [0.0, 0.0, 0.0];
    }
    else {
        var I_r = [0.5 / (lightDis * lightDis), 0.5 / (lightDis * lightDis), 0.5 / (lightDis * lightDis)];
    }
    // console.log("I_r: ", I_r);
    var diffuseFactor = Math.max(dot(v_normal, lightVec), 0.0);
    var diffuseReflection = vecscalarMul(vecMul(kd, I_r), m);
    // console.log("diffuseReflection: ", diffuseReflection)
    diffuseReflection = vecscalarMul(diffuseReflection, diffuseFactor * mass * 100.0);

    //Specular Highlights
    var store = vecscalarMul(I_r, ks);
    var v = normalize(vecSub(camPos, loc));
    var h = normalize(vecAdd([lightVec, v]));

    var max_val = Math.pow(Math.max(0.0, dot(v_normal_normalize, h)), p);
    var specHighlight = vecscalarMul(store, max_val);
    specHighlight = vecscalarMul(specHighlight, mass * 1.0 * norm(velocity));

    // console.log("ambLighting: ", ambLighting);
    // console.log("diffuseReflection: ", diffuseReflection);
    // console.log("specHighlight: ", specHighlight);
    outColor = vecAdd([ambLighting, diffuseReflection, specHighlight]); //, diffuseReflection, specHighlight
    outColor = [outColor[0] * 255, outColor[1] * 255, outColor[2] * 255, 255]
    return outColor;
}

function draw() {
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

    // heat --
    minT = f.temp[0];
    maxT = f.temp[0];
    for (var i = 0; i < f.numCells; i++) {
        minT = Math.min(minT, f.temp[i]);
        maxT = Math.max(maxT, f.temp[i]);
    }
    // -- heat

    id = c.getImageData(0, 0, canvas.width, canvas.height)

    var color = [255, 255, 255, 255]
    for (var i = 0; i < f.numX; i++) {
        for (var j = 0; j < f.numY; j++) {
            if (scene.showPressure) {
                var p = f.p[i * n + j];
                var s = f.m[i * n + j];
                color = getSciColor(p, minP, maxP);
                if (scene.showSmoke) {
                    color[0] = Math.max(0.0, color[0] - 255 * s);
                    color[1] = Math.max(0.0, color[1] - 255 * s);
                    color[2] = Math.max(0.0, color[2] - 255 * s);
                }
            }
            else if (scene.showTemp) {
                var t = f.temp[i * n + j];
                color = getSciColor(t, minT, maxT);
                if (scene.showSmoke) {
                    color[0] = Math.max(0.0, color[0]);
                    color[1] = Math.max(0.0, color[1]);
                    color[2] = Math.max(0.0, color[2]);
                }
            }
            else if (scene.showSmoke) {
                if (scene.smokeAppOne) {
                    var s = f.m[i * n + j];
                    color[0] = 255 * s;
                    color[1] = 255 * s;
                    color[2] = 255 * s;
                    color[3] = 255;

                } else if (scene.smokeAppTwo) {
                    mass = f.m[i * n + j]
                    velocity = [f.u[i * n + j], f.v[i * n + j]]
                    loc = [i, j, 0.0]
                    var colorValue = phongShader(mass, velocity, loc);
                    color[0] = colorValue[0];
                    color[1] = colorValue[1];
                    color[2] = colorValue[2];
                    color[3] = colorValue[3];

                } else if (scene.smokeAppThree) {
                    mass = f.m[i * n + j];
                    var colorValue = vangoghShader(mass);
                    color[0] = colorValue[0];
                    color[1] = colorValue[1];
                    color[2] = colorValue[2];
                    color[3] = colorValue[3];

                } else if (scene.fire) {
                    var heat = f.temp[i * n + j];
                    var colorValue = fireShader(heat);
                    color[0] = colorValue[0];
                    color[1] = colorValue[1];
                    color[2] = colorValue[2];
                    color[3] = colorValue[3];
                } else if (scene.blocky) {
                    var mass = f.m[i * n + j];
                    var colorValue = blockyShader(mass, [0, 255, 255, 255]);
                    color[0] = colorValue[0];
                    color[1] = colorValue[1];
                    color[2] = colorValue[2];
                    color[3] = colorValue[3];
                } else if (scene.color) {
                    var mass = f.m[i * n + j];
                    var colorValue = colorShader(mass, [255, 255, 0, 255]);
                    color[0] = colorValue[0];
                    color[1] = colorValue[1];
                    color[2] = colorValue[2];
                    color[3] = colorValue[3];
                } else if (scene.realistic) {
                    var colorValue = realisticShader(c.getImageData(0, 0, canvas.width, canvas.height), f);
                    color[0] = colorValue[0];
                    color[1] = colorValue[1];
                    color[2] = colorValue[2];
                    color[3] = colorValue[3];
                }

                if (scene.sceneNr == 2) {
                    color = getSciColor(s, 0.0, 1.0);
                }
            }
            else if (f.s[i * n + j] == 0.0) {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }

            var x = Math.floor(cX(i * h));
            var y = Math.floor(cY((j + 1) * h));
            var cx = Math.floor(cScale * cellScale * h) + 1;
            var cy = Math.floor(cScale * cellScale * h) + 1;

            r = color[0];
            g = color[1];
            b = color[2];
            a = color[3]; //added

            for (var yi = y; yi < y + cy; yi++) {
                var p = 4 * (yi * canvas.width + x)

                for (var xi = 0; xi < cx; xi++) {
                    id.data[p++] = r;
                    id.data[p++] = g;
                    id.data[p++] = b;
                    // id.data[p++] = 255;
                    id.data[p++] = a;
                }
            }
        }
    }

    c.putImageData(id, 0, 0);

    if (scene.showVelocities) {

        c.strokeStyle = "#000000";
        scale = 0.02;

        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {

                var u = f.u[i * n + j];
                var v = f.v[i * n + j];

                c.beginPath();

                x0 = cX(i * h);
                x1 = cX(i * h + u * scale);
                y = cY((j + 0.5) * h);

                c.moveTo(x0, y);
                c.lineTo(x1, y);
                c.stroke();

                x = cX((i + 0.5) * h);
                y0 = cY(j * h);
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
                    l = Math.sqrt(u * u + v * v);
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
        else if (scene.showTemp)
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

    if (scene.showGrenade) {
        c.strokeW
        r = scene.grenadeRadius;
    }

    if (scene.showPressure) {
        var s = "pressure: " + minP.toFixed(0) + " - " + maxP.toFixed(0) + " N/m";
        c.fillStyle = "#000000";
        c.font = "16px Arial";
        c.fillText(s, 10, 35);
    }

    if (scene.showTemp) {
        var s = "temperature: " + minT.toFixed(0) + " - " + maxT.toFixed(0) + " ";
        c.fillStyle = "#000000";
        c.font = "16px Arial";
        c.fillText(t, 10, 35);
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

    for (var i = 1; i < f.numX - 2; i++) {
        for (var j = 1; j < f.numY - 2; j++) {

            f.s[i * n + j] = 1.0;

            dx = (i + 0.5) * f.h - x;
            dy = (j + 0.5) * f.h - y;

            if (dx * dx + dy * dy < r * r) {
                f.s[i * n + j] = 0.0;
                if (scene.sceneNr == 2)
                    f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                else
                    f.m[i * n + j] = 1.0;
                f.u[i * n + j] = vx;
                f.u[(i + 1) * n + j] = vx;
                f.v[i * n + j] = vy;
                f.v[i * n + j + 1] = vy;
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


    for (var i = 1; i < f.numX - 2; i++) {
        for (var j = 1; j < f.numY - 2; j++) {
            if (true) {
                f.s[i * n + j] = 1.0;


                dx = (i + 0.5) * f.h - x;
                dy = (j + 0.5) * f.h - y;


                if (dx * dx + dy * dy < r * r) {
                    f.s[i * n + j] = 0.0;
                    if (scene.sceneNr == 2)
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                    else
                        if (f.m[i * n + j] < 0.1) {
                            f.m[i * n + j] = 0.1;
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


                    f.u[(i - 1) * n + (j - 1)] += strength;
                    f.v[(i - 1) * n + (j - 1)] += strength;


                    f.u[(i - 1) * n + (j)] += strength;
                    f.v[(i - 1) * n + (j)] += strength;


                    f.u[(i - 1) * n + (j + 1)] += strength;
                    f.v[(i - 1) * n + (j + 1)] += strength;


                    f.u[(i) * n + (j - 1)] += strength;
                    f.v[(i) * n + (j - 1)] += strength;


                    // f.u[(i)*n + (j)] += 10;
                    // f.v[(i)*n + (j)] += 10;


                    f.u[(i + 1) * n + (j + 1)] += strength;
                    f.v[(i + 1) * n + (j + 1)] += strength;


                    f.u[(i + 1) * n + (j - 1)] += strength;
                    f.v[(i + 1) * n + (j - 1)] += strength;


                    f.u[(i + 1) * n + (j)] += strength;
                    f.v[(i + 1) * n + (j)] += strength;


                    f.u[(i + 1) * n + (j + 1)] += strength;
                    f.v[(i + 1) * n + (j + 1)] += strength;
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
        for (var i = 1; i < f.numX - 2; i++) {
            for (var j = 1; j < f.numY - 2; j++) {
                f.s[i * n + j] = 1.0;

                dx = (i + 0.5) * f.h - flare_x;
                dy = (j + 0.5) * f.h - flare_y;

                if (dx * dx + dy * dy < r * r) {
                    f.s[i * n + j] = 0.0;
                    if (scene.sceneNr == 2)
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                    else
                        if (f.m[i * n + j] > 0.2) {
                            f.m[i * n + j] = 0.5;
                        }
                    f.u[i * n + j] = 1 / 25;
                    f.u[(i + 1) * n + j] = 3 / 25;
                    f.v[i * n + j] = 5 / 25;
                    f.v[i * n + j + 1] = 7 / 25;
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

    for (var i = 1; i < f.numX - 2; i++) {
        for (var j = 1; j < f.numY - 2; j++) {
            // console.log("f.m[i*n + j]: ", f.m[i*n +j])
            if (true) {
                f.s[i * n + j] = 1.0;

                dx = (i + 0.5) * f.h - x;
                dy = (j + 0.5) * f.h - y;


                if (dx * dx + dy * dy < r * r) {
                    f.s[i * n + j] = 0.0;
                    if (scene.sceneNr == 2)
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr)
                    else
                        if (f.m[i * n + j] > 0.2) {
                            f.m[i * n + j] = 0.85;
                        }
                    f.u[i * n + j] = 0.005;
                    f.u[(i + 1) * n + j] = 0.005;
                    f.v[i * n + j] = 0.005;
                    f.v[i * n + j + 1] = 0.005;
                }
                else if (dx * dx + dy * dy < r * r + (0.00325) && dx * dx + dy * dy > r * r - (0.00325)) {
                    f.m[i * n + j] = 1.0;
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

    for (var i = 1; i < f.numX - 2; i++) {
        for (var j = 1; j < f.numY - 2; j++) {
            f.s[i * n + j] = 1.0;

            var clear = true;

            var num_activeSmokes = scene.activeSmokes.length;
            if (num_activeSmokes == 1) {
                f.m[i * n + j] = 1.0;
            }
            for (var k = 0; k < num_activeSmokes; k++) {
                var smoke_x = scene.activeSmokes[k][0];
                var smoke_y = scene.activeSmokes[k][1];

                dx = (i + 0.5) * f.h - smoke_x;
                dy = (j + 0.5) * f.h - smoke_y;

                var store = dx * dx + dy * dy > r * r + (0.004);
                clear = clear && store;
            }

            if (clear) {
                f.m[i * n + j] = 1.0;
            }
        }
    }
}

function tempReductionOverTime() {
    var n = this.numY;
    var reductionConst = 0.02;
    for (var i = 1; i < this.numX - 1; i++) {
        for (var j = 1; j < this.numY - 1; j++) {
            this.temp[i * n + j] += (0.3 - this.temp[i(n + j)]) * reductionConst;
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

    setObstacle(x, y, true);
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
            setObstacle(x, y, false);
        }
        // if (scene.flaregrenade) {
        // 	setFlareGrenade(x, y);
        // }
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
    else if (scene.flaregrenade) {
        startFlareGrenade(event.x, event.y);
        console.log("scene.activeFlares: ", scene.activeFlares);
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
}, { passive: false });


document.addEventListener('keydown', event => {
    switch (event.key) {
        case 'p': scene.paused = !scene.paused; break;
        case 'm': scene.paused = false; simulate(); scene.paused = true; break;
    }
});

function toggleStart() {
    var button = document.getElementById('startButton');
    if (scene.paused)
        button.innerHTML = "Stop";
    else
        button.innerHTML = "Start";
    scene.paused = !scene.paused;
}

// main -------------------------------------------------------

function simulate() {
    if (!scene.paused) {
        scene.time += scene.dt;
        scene.fluid.simulate(scene.dt, scene.gravity, scene.numIters);
        scene.frameNr++;
    }

    if (scene.activeFlares.length > 0) {
        displayFlareGrenades();
    }

    if (!scene.ball) {
        setObstacle(0.4, 0.5, true)
    }
    if (scene.activeSmokes.length > 0) {
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