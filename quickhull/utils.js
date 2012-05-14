function getRandomPoints3D(numPoint, xMax, yMax, zMax) {
    var points = new Array();
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ xMax /4 + r * Math.cos(theta), yMax/2 + 2 * r * Math.sin(theta),  zMax/5 + 2 * r * Math.sin(theta)] );
    }
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ xMax /4 * 3 +  r * Math.cos(theta), yMax/2 +  r * Math.sin(theta),  zMax /5 * 3 +  r * Math.cos(theta) ] );
    }
    return points;
}

function getRandomPoints2D(numPoint, xMax, yMax) {
    var points = new Array();
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ (xMax /4 + r * Math.cos(theta)) * Math.pow(-1, Math.round(Math.random())), 
                       (yMax/2 + 2 * r * Math.sin(theta)) * Math.pow(-1, Math.round(Math.random()))] );
    }
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ (xMax /4 * 3 +  r * Math.cos(theta))* Math.pow(-1, Math.round(Math.random())), 
                       (yMax/2 +  r * Math.sin(theta))* Math.pow(-1, Math.round(Math.random()))] );
    }
    return points;
}

function random2D(numPoints, maxX, maxY) {
	var points = new Array();

	for ( var i = 0; i < numPoints; i++) {
		
		points.push([Math.floor(Math.random() * maxX) * Math.pow(-1, Math.round(Math.random())), 
		             Math.floor(Math.random() * maxY) * Math.pow(-1, Math.round(Math.random()))]);
	}
	
	return points;
}

function random3D(numPoints, maxX, maxY, maxZ) {
	var points = new Array();

	for ( var i = 0; i < numPoints; i++) {
		
		points.push([Math.random() * maxX * Math.pow(-1, Math.round(Math.random())), 
		             Math.random() * maxY * Math.pow(-1, Math.round(Math.random())),
		             Math.random() * maxZ * Math.pow(-1, Math.round(Math.random()))]);
	}
	
	return points;
}

function random4D(numPoints, maxX, maxY, maxZ, maxW) {
	var points = new Array();

	for ( var i = 0; i < numPoints; i++) {
		
		points.push([Math.random() * maxX * Math.pow(-1, Math.round(Math.random())), 
		             Math.random() * maxY * Math.pow(-1, Math.round(Math.random())),
		             Math.random() * maxZ * Math.pow(-1, Math.round(Math.random())),
		             Math.random() * maxW * Math.pow(-1, Math.round(Math.random()))]);
	}
	
	return points;
}

function random4DPos(numPoints, maxX, maxY, maxZ, maxW) {
	var points = new Array();

	for ( var i = 0; i < numPoints; i++) {
		
		points.push([Math.floor(Math.random() * maxX) + 1, 
		             Math.floor(Math.random() * maxY) + 1,
		             Math.floor(Math.random() * maxZ) + 1,
		             Math.floor(Math.random() * maxW) + 1]);
	}
	
	return points;
}

function convHullToPlasm(vertices, facets) {
	var tmpVertices;
	var cells = new Array();
	for ( var i = 0; i < facets.length; i++) {
		tmpVertices = facets[i].vertices;
		cells.push([]);
		for ( var j = 0; j < tmpVertices.length; j++) {
			for ( var k = 0; k < vertices.length; k++) {
				if ((quickhull._utils.comparePoints(vertices[k], tmpVertices[j]) === 0)) {
					cells[i].push(k);
				}
			}
		}
	}

	return cells;
}

var toCells = function(vertices, facets) {
	var cells = new Array();
	facets.forEach(function(facet, i){
		cells.push([]);
		vertices.forEach(function(vertex,j){
			facet.vertices.forEach(function(facetVert,k){
				if ((quickhull._utils.comparePoints(facetVert, vertex) === 0)) {
					cells[i].push(j);
				}
			});
		});
	});
	return cells;
};

var schlegel3D = function (d) {
	return function(point){
		var denominator = point[3]/d;
		return [point[0]/denominator, point[1]/denominator, point[2]/denominator];
	};
};

function qhPlotPoints(pts) {
    ctx = document.getElementById('qh_demo').getContext('2d');
    ctx.clearRect(0,0,400,400);
    ctx.fillStyle = 'rgb(0,0,0)';
    for (var idx in pts) {
        var pt = pts[idx];
        ctx.fillRect(pt[0],pt[1],2,2);
    }
}

function qhPlotPoint(pt) {
    ctx = document.getElementById('qh_demo').getContext('2d');
    ctx.fillStyle = 'rgb(255,0,0)';
    
    ctx.fillRect(pt[0],pt[1],2,2);
    
}

function plotBaseLine(baseLine,color) {
    var ctx = document.getElementById('qh_demo').getContext('2d');
    var pt1 = baseLine[0];
    var pt2 = baseLine[1];
    ctx.save();
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(pt1[0],pt1[1]);
    ctx.lineTo(pt2[0],pt2[1]);
    ctx.closePath();
    ctx.stroke();
    ctx.restore();
} 

Array.prototype.shuffle = function() {
	for(var j, x, i = this.length; i; j = Math.floor(Math.random() * i), x = this[--i], this[i] = this[j], this[j] = x);
};

Array.prototype.mix = function() {
	var ret = new Array();
	while (this[0]) {
		ret.push(
			this.splice(
				parseInt(Math.random()*this.length),
				1)[0]
			);
	}
	while (ret[0]) {
		this.push(ret.shift());
	}
	return this;
};