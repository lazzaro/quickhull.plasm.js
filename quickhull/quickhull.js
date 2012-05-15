/*!
 * quickhull.js
 * Quickhull algoritm in javascript translated from 
 * http://www.qhull.org
 * 
 */
!(function (exports) {

	/**
	 * Variables.
	 */

	var abs = Math.abs;
	var sqrt = Math.sqrt;
	
	var facetId = 0; // id globale per le facce
	var idVisit = 0; // id di visita per identificare la singola iterazione
	
	var displayGlob = 0;
	var logGlob = 0;
	var steps;
	
	// Oggetto per la memorizzazione della struttura della singola faccia
	function Facet(){
		this.id = facetId++;
		this.neighbors = new Array();
		this.furthestDist = 0.0;
		this.furthestPoint = null;
		this.vertices = new Array();
		this.topoOriented = true;
		this.normal = new Array();
		this.offset = 0.0;
		this.outsideSet = new Array();
		this.visitId = idVisit;
		this.visible = false;
	};

	/**
	 * Library namespace.
	 */

	var quickhull = exports.quickhull = {};

	/**
	 * Library version.
	 */

	quickhull.version = '0.2.0';

	/**
	 * utils namespace
	 * @api private
	 */

	quickhull._utils = {};



	/**
	 * verify if a point is already inside the set
	 * 
	 * @param {Array|Float32Array} set - set of points
	 * @param {Array|Float32} point - point to verify
	 * @return {Boolean} - if the point is already in
	 * 
	 */

	var isIn = quickhull._utils.isIn = function(set, point) {
		var isIn = false, equal;
		var tempPoint;

		for ( var i = 0; i < set.length && !isIn; i++) {
			tempPoint = set[i];
			
			equal = (comparePoints(tempPoint, point) === 0);

			if (equal) {
				isIn = true;
			}
		}

		return isIn;
	};

	/**
	 * Add a point to a set only if not already present
	 * 
	 * @param {Array|Float32Array} set - set of points
	 * @param {Array|Float32} point - point to add
	 * 
	 */

	var addNoDup = quickhull._utils.addNoDup = function (set, point){
		
		if (!isIn(set, point)) {
			set.push(point);
		}
	};
	
	/**
	 * Remove a point from a set - the set not contain duplicate points
	 * 
	 * 
	 * @param {Array|Float32Array} set - set of points
	 * @param {Array|Float32} point - point to remove
	 * 
	 */

	var removePoint = quickhull._utils.removePoint = function (set, point){
		
		for ( var i = 0; i < set.length; i++) {
			if ((comparePoints(set[i], point) === 0)) {
				set.splice(i,1);
			}
		}
	};
	
	/**
	 * Determinant of 2x2 matrix
	 * 
	 * @param {Number} a1 - coordinate 0,0
	 * @param {Number} a2 - coordinate 0,1
	 * @param {Number} b1 - coordinate 1,0
	 * @param {Number} b2 - coordinate 1,1
	 * @return {Number} determinant
	 * 
	 */
	
	var det2 = quickhull._utils.det2 = function(a1, a2, b1, b2) {
		return (a1 * b2) - (a2 * b1);
	};
	
	/**
	 * Determinant of 3x3 matrix
	 * 
	 * @param {Number} a1 - coordinate 0,0
	 * @param {Number} a2 - coordinate 0,1
	 * @param {Number} a2 - coordinate 0,2
	 * @param {Number} b1 - coordinate 1,0
	 * @param {Number} b2 - coordinate 1,1
	 * @param {Number} b2 - coordinate 1,2
	 * @param {Number} c1 - coordinate 2,0
	 * @param {Number} c2 - coordinate 2,1
	 * @param {Number} c2 - coordinate 2,2
	 * @return {Number} determinant
	 * 
	 */
	
	var det3 = quickhull._utils.det3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3) {
		return (a1) * det2(b2,b3,c2,c3) - (b1) * det2(a2,a3,c2,c3) + (c1)*det2(a2,a3,b2,b3);
	};

	/**
	 * Compute the determinant of a simplex
	 * 
	 * @param {Array|Float32} point - point apex
	 * @param {Array|Float32Array} simplex - base points
	 * @param {Number} dim - dimension
	 */

	var detSimplex = quickhull._utils.detSimplex = function(point, simplex, dim) {
		
		var j = 0, det;
		var sign;
		var rows = new Array(), rowsGauss = new Array(), coordPointSimpl, row = new Array();
		
		for ( var i = 0; i < simplex.length && !(j === dim); i++) {
			coordPointSimpl = simplex[i];
			for ( var k = 0; k < dim; k++) {
				row.push(coordPointSimpl[k] - point[k]);
			}
			
			rows.push(row);
			row = new Array();
			j = rows.length;
		}
		
		
		// calcolo del determinante
		//TODO implementare gestione nearzero
		if (dim === 2) {
			det = det2(rows[0][0], rows[0][1], rows[1][0], rows[1][1]); 
		} else if (dim === 3) {
			det = det3(rows[0][0], rows[0][1], rows[0][2], rows[1][0], rows[1][1], rows[1][2], rows[2][0], rows[2][1], rows[2][2]);
		} else {
			sign = gaussElim(rows, false);
			
			det = 1.0;
			for (var l = 0; l < 0; l++) {
				det = det * rowsGauss[l][l];
			}
			
			if (sign) {
				det = -det;
			}
		}
		
		return det;
	};
	
	/**
	 * Gaussian elimination with partial pivoting (perform gaussian elimination step)
	 * 
	 * @param {Array|Float32Array} rows - the matrix's rows
	 * @return {Boolean} sign - if the determinant must change sign
	 * 
	 */
	var gaussElim = quickhull._utils.gaussElim = function(rows, sign) {
		var numRows = rows.length;
		var numCols = rows[0].length;
		var pivotAbs, pivoti, temp, rowP, pivot;
		var n;
		
		for ( var k = 0; k < numRows; k++) {
			pivotAbs = abs(rows[k][k]); 
			pivoti = k;
			
			for ( var i = k + 1; i < numRows; i++) {
				temp = abs(rows[i][k]);
				if (temp > pivotAbs) {
					pivotAbs = temp;
					pivoti = i;
				}
			}
			
			if (pivoti != k) {
				rowP = rows[pivoti];
				rows[pivoti] = rows[k];
				rows[k] = rowP;
				sign = !sign;
			}
			
			pivot = rows[k][k];
			//TODO gestione nearzero
			if(pivot !== 0){				
				for ( var i = k + 1; i < numRows; i++) {
					n = rows[i][k] / pivot;
					
					for ( var j = k - 1; j < numCols; j++) {
						
						rows[i][j] = rows[i][j] - n * rows[k][j];
					}
				}
			}
		}
		
		return sign;
	};
	
	/**
	 * constructs the initial hull
	 * 
	 * @param {Array|Float32Array} vertices - initial vertices
	 * @return {Array|Facet} initial hull
	 * 
	 */

	var initialHull = quickhull._utils.initialHull = function (vertices) {
		var newFacets = new Array();
		var newFacet;
		var interiorPoint;
		var dist;
		var cells, cell, conv;
		
		if (displayGlob === 2) {
			cells = new Array();
		}
		
		interiorPoint = getCenter(vertices);
		
		for ( var i = 0; i < vertices.length; i++) {
			newFacet = new Facet();
			
			if (displayGlob === 2) {
				cell = new Array();
			}
			
			for ( var j = 0; j < vertices.length; j++) {
				var k = (vertices.length - 1 - i);
				if (j != k) {
					newFacet.vertices.push(vertices[j]);
					
					if (displayGlob === 2) {
						cell.push(j);
					}
				}
			}
			
			if (displayGlob === 2) {
				cells.push(cell);
			}
			
			newFacets.push(newFacet);
			
			setFacetPlane(newFacet);
			
			dist = distPlane(interiorPoint, newFacet);
			
			if (dist > 0) {
				for ( var k = 0; k < newFacet.normal.length; k++) {
					newFacet.normal[k] = -newFacet.normal[k];
				}
				newFacet.offset = -newFacet.offset;
				newFacet.topoOriented = false;
			}
		}
		
		for ( var i = 0; i < newFacets.length; i++) {
			newFacet = newFacets[i];
			for ( var j = 0; j < newFacets.length; j++) {
				if (newFacet.id !== newFacets[j].id) {
					newFacet.neighbors.push(newFacets[j]);
				}
			}
		}
		//TODO implementare checkflipped, gestione minangle, qh_checkpolygon, qh_checkconvex
		
		return newFacets;
	};
	
	/**
	 * return distance from point to facet
	 * 
	 * @param {Array|Float32Array} point - point
	 * @param {Facet} facet - facet
	 * @return {Number} distance of point from facet
	 */

	var distPlane = quickhull._utils.distPlane = function (point, facet) {
		var dist = facet.offset, normal = facet.normal;
		
		for ( var i = 0; i < normal.length; i++) {
			dist = dist + (point[i]*normal[i]);
		}
		
		return dist;
	};
	
	/**
	 * returns arithmetic center of a set of vertices
	 * 
	 * @param {Array|Float32Array} vertices - set of vertices
	 * @return {Float32Array} center point
	 */

	var getCenter = quickhull._utils.getCenter = function (vertices) {
		
		var center = new Array();
		var size = vertices.length;
		var dim = vertices[0].length;
		var coordinata;
		
		//TODO implementare controllo che siano almeno due punti
		for ( var i = 0; i < dim; i++) {
			coordinata = 0.0;
			for ( var j = 0; j < vertices.length; j++) {
				coordinata = coordinata + vertices[j][i];
			}
			center.push(coordinata/size);
		}
		
		return center;
	};
	
	/**
	 * sets the hyperplane for a facet and the offset
	 * 
	 * @param {Facet} facet - the facet to set hyperplane
	 * 
	 */

	var setFacetPlane = quickhull._utils.setFacetPlane = function (facet) {
		var point0 = facet.vertices[0], normal = new Array();
		var offset, rows = new Array();
		
		if (point0.length <= 4) {
			for ( var i = 0; i < facet.vertices.length; i++) {
				rows.push(facet.vertices[i]);
			}
			
			normal = hyperplaneDet(rows, point0, facet.topoOriented);
		}
		
		if (point0.lenght > 4) { //TODO || nearzero
			rows = new Array();
			
			for ( var i = 0; i < facet.vertices.length; i++) {
				rows.push(facet.vertices[i]);
			}
			
			normal = hyperplaneGauss(rows, point0, facet.topoOriented);
		}
		
		facet.normal = normal;
		offset = -(point0[0] * normal[0]);
		
		for ( var i = 1; i < normal.length; i++) {
			offset = offset - (point0[i] * normal[i]);
		}
		
		facet.offset = offset;
	};
	
	/**
	 * given an upper-triangular rows array and a sign,
     * solve for normal equation x using back substitution over rows U
	 * 
	 * @param {Array|Float32Array} rows - one row per point
	 * @param {Number} numrow - number of rows
	 * @param {Number} numcol - number of cols
	 * @param {Boolean} sign - 
	 * @return {Float32Array} normal normalized
	 * 
	 */

	var backNormal = quickhull._utils.backNormal = function (rows, numrow, numcol, sign) {
		var normal = new Array(numcol), diagonal;
		normal[normal.length - 1] = (sign ? -1.0 : 1.0);
		for ( var i = numrow - 2; i >= 0; i--) {
			normal[i] = 0.0;
			for ( var j = i + 1; j < numcol; j++) {
				normal[i] = normal[i] - (rows[i][j] * normal[j]);
			}
			diagonal = rows[i][i];
			//TODO gestire caso diagonal uguale a zero
			normal[i] = normal[i]/diagonal;
		}
		
		return normal;
	};
	
	/**
	 * return normalized hyperplane equation from oriented simplex
	 * 
	 * @param {Array|Float32Array} rows - one row per point
	 * @param {Float32Array} point0 - any row
	 * @param {Boolean} topoOrient - flips all signs
	 * @return {Float32Array} normal normalized
	 * 
	 */

	var hyperplaneGauss = quickhull._utils.hyperplaneGauss = function (rows, point0, topoOrient) {
		var dim = rows[0].length;
		var normal;
		var sign = gaussElim(rows, topoOrient);
		
		for ( var i = 0; i < dim; i++) {
			if (rows[i][i] < 0) {
				sign = !sign;
			}
		}
		
		normal = backNormal(rows, rows.length - 1, rows[0].length, sign);
		//TODO gestione nearzero
		
		normalize2(normal, dim, topoOrient);
		
		return normal;
	};
	
	/**
	 * return normalized hyperplane equation from oriented simplex
	 * 
	 * @param {Array|Float32Array} rows - one row per point
	 * @param {Float32Array} point0 - any row
	 * @param {Boolean} topoOrient - flips all signs
	 * @return {Float32Array} normal normalized
	 * 
	 */

	var hyperplaneDet = quickhull._utils.hyperplaneDet = function (rows, point0, topoOrient) {
		/*given two indices into rows[],

         compute the difference between X, Y, or Z coordinates
		  #define dX( p1,p2 )  ( *( rows[p1] ) - *( rows[p2] ))
		  #define dY( p1,p2 )  ( *( rows[p1]+1 ) - *( rows[p2]+1 ))
		  #define dZ( p1,p2 )  ( *( rows[p1]+2 ) - *( rows[p2]+2 ))
		  #define dW( p1,p2 )  ( *( rows[p1]+3 ) - *( rows[p2]+3 ))*/
		var dim = rows[0].length;
		var normal = new Array();
		
		if (dim === 2) {
			normal.push(rows[1][1] - rows[0][1]);//dY(1,0)
			normal.push(rows[0][0] - rows[1][0]);//dX(0,1)
			normalize2(normal, dim, topoOrient);
		} else if (dim === 3) {
			normal.push(det2(rows[2][1] - rows[0][1],//dY(2,0) 
					         rows[2][2] - rows[0][2],//dZ(2,0)
					         rows[1][1] - rows[0][1],//dY(1,0)
					         rows[1][2] - rows[0][2]));//dZ(1,0)
			
			normal.push(det2(rows[1][0] - rows[0][0],//dX(1,0)
			         		 rows[1][2] - rows[0][2],//dZ(1,0)
			         		 rows[2][0] - rows[0][0],//dX(2,0)
			         		 rows[2][2] - rows[0][2]));//dZ(2,0)
			
			normal.push(det2(rows[2][0] - rows[0][0],//dX(2,0)
			         		 rows[2][1] - rows[0][1],//dY(2,0)
			         		 rows[1][0] - rows[0][0],//dX(1,0)
			         		 rows[1][1] - rows[0][1]));//dY(1,0)
			
			normalize2(normal, dim, topoOrient);
		} else if (dim === 4) {
			normal.push(-det3(rows[2][1] - rows[0][1],//dY(2,0)
	         		 		  rows[2][2] - rows[0][2],//dZ(2,0)
	         		 		  rows[2][3] - rows[0][3],//dW(2,0)
	         		 		  rows[1][1] - rows[0][1],//dY(1,0)
	         		 		  rows[1][2] - rows[0][2],//dZ(1,0)
	         		 		  rows[1][3] - rows[0][3],//dW(1,0)
	         		 		  rows[3][1] - rows[0][1],//dY(3,0)
	         		 		  rows[3][2] - rows[0][2],//dZ(3,0)
	         		 		  rows[3][3] - rows[0][3]));//dW(3,0)
			
			normal.push(det3(rows[2][0] - rows[0][0],//dX(2,0) 
   		 		  			 rows[2][2] - rows[0][2],//dZ(2,0)
   		 		  			 rows[2][3] - rows[0][3],//dW(2,0)
   		 		  			 rows[1][0] - rows[0][0],//dX(1,0) 
   		 		  			 rows[1][2] - rows[0][2],//dZ(1,0)
   		 		  			 rows[1][3] - rows[0][3],//dW(1,0)
   		 		  			 rows[3][0] - rows[0][0],//dX(3,0) 
   		 		  			 rows[3][2] - rows[0][2],//dZ(3,0)
   		 		  			 rows[3][3] - rows[0][3]));//dW(3,0)
			
			normal.push(-det3(rows[2][0] - rows[0][0],//dX(2,0) 
   		 		  			  rows[2][1] - rows[0][1],//dY(2,0)
   		 		  			  rows[2][3] - rows[0][3],//dW(2,0)
   		 		  			  rows[1][0] - rows[0][0],//dX(1,0) 
   		 		  			  rows[1][1] - rows[0][1],//dY(1,0)
   		 		  			  rows[1][3] - rows[0][3],//dW(1,0)
   		 		  			  rows[3][0] - rows[0][0],//dX(3,0) 
   		 		  			  rows[3][1] - rows[0][1],//dY(3,0)
   		 		  			  rows[3][3] - rows[0][3]));//dW(3,0)
			
			normal.push(det3(rows[2][0] - rows[0][0],//dX(2,0)
	 		  			  	 rows[2][1] - rows[0][1],//dY(2,0)
	 		  			  	 rows[2][2] - rows[0][2],//dZ(2,0)
	 		  			  	 rows[1][0] - rows[0][0],//dX(1,0)
	 		  			  	 rows[1][1] - rows[0][1],//dY(1,0)
	 		  			  	 rows[1][2] - rows[0][2],//dZ(1,0)
	 		  			  	 rows[3][0] - rows[0][0],//dX(3,0)
	 		  			  	 rows[3][1] - rows[0][1],//dY(3,0)
	 		  			  	 rows[3][2] - rows[0][2]));//dZ(3,0)
			
			normalize2(normal, dim, topoOrient);
		}
		return normal;
	};
	
	/**
	 * normalize a vector
	 * 
	 * @param {Float32Array} normal - vector to normalize
	 * @param {Boolean} topoOrient - flips all signs
	 * 
	 */
	
	var normalize2 = quickhull._utils.normalize2 = function (normal, dim, topoOrient) {
		var somma = 0.0, norm = 0.0, temp;
		for ( var i = 0; i < dim; i++) {
			somma = somma + (normal[i] * normal[i]);
		}
		
		norm = sqrt(somma);
		
		if (norm === 0.0) {
			temp = sqrt(1.0/dim);
			for ( var i = 0; i < dim; i++) {
				normal[i] = temp;
			}
		} else {
			if (!topoOrient) {
				norm = -norm;
			}
			for ( var i = 0; i < dim; i++) {
				//TODO gestire divisioni per zero
				normal[i] = normal[i] / norm;
			}
		}
	};
	
	/**
	 * Compare points with same dimension
	 * 
	 * @param {Float32Array} point1 - first point
	 * @param {Float32Array} point2 - second point
	 * @return {Boolean} if equal
	 * 
	 */ 

	var comparePoints = quickhull._utils.comparePoints = function(point1, point2) {
		var compare = 0;
		for ( var i = 0; i < point1.length && (compare === 0); i++) {
			if (point1[i] > point2[i]) {
				compare = 1;
			} else if (point1[i] < point2[i]) {
				compare = -1;
			}
		}
		
		return compare;
	};
	
	/**
	 * if outside add point to facet's outside set
	 * 
	 * @param {Array|Float32Array} point - the point to add if outside of facet
	 * @param {Facet} facet - the facet
	 * @return {Number} distance of point from facet
	 * 
	 */
	
	var addOutside = quickhull._utils.addOutside = function(point, facet) {
		var dist;
		
		dist = distPlane(point, facet);
		
		if (dist > 0) {
			if (facet.furthestPoint === null) {
				facet.furthestPoint = point;
				facet.furthestDist = dist;
			} else if (dist > facet.furthestDist) {
				facet.outsideSet.push(facet.furthestPoint);
				facet.furthestPoint = point;
				facet.furthestDist = dist;
			} else {
				facet.outsideSet.push(point);
			}
		}
		return dist;
	};
	
	/**
	 * partitions all points in points/numpoints to the outsidesets of facets
	 * 
	 * @param {Array|Facet} facets - set of facets
	 * @param {Array|Float32Array} vertices - set of points facets' vertices
	 * @param {Array|Float32Array} points - set of points
	 * @param {Number} - number of points
	 * 
	 */ 
	
	var partitionAll = quickhull._utils.partitionAll = function(facets, vertices, points) {
		for ( var i = 0; i < facets.length; i++) {
			for ( var j = 0; j < points.length; j++) {
				if (!isIn(vertices, points[j])) {
					addOutside(points[j], facets[i]);
				}
			}
		}
	};
	
	/**
	 * return the facet with furthest of all furthest points searches all facets
	 * 
	 * @param {Array|Facet} facets - set of facets
	 * @return {Facet} the facet with furthest 
	 * 
	 */ 
	
	var furthestNext = quickhull._utils.furthestNext = function(facets) {
		var bestFacet = facets[0], bestDist = bestFacet.furthestDist;
		
		for ( var i = 1; i < facets.length; i++) {
			if (bestDist < facets[i].furthestDist) {
				bestDist = facets[i].furthestDist;
				bestFacet = facets[i];
			}
		}
		
		if (bestDist === 0.0) {
			bestFacet = null;
		}
		
		return bestFacet;
	};
	
	/**
	 * given a visible facet, find the point's horizon and visible facets for all facets
	 * 
	 * @param {Array} point - the point to find horizon
	 * @param {Facet} facet - the first facet
	 * @return {Array|Facet} the visible facets
	 * 
	 */ 
	
	var findHorizon = quickhull._utils.findHorizon = function(point, facet) {
		idVisit++;
		var visible = new Array();
		var tmpNeighbors, neighbor, dist;
		
		facet.visitId = idVisit;
		facet.visible = true;
		
		visible.push(facet);
		
		for ( var i = 0; i < visible.length; i++) {
			tmpNeighbors = visible[i].neighbors;
			
			for ( var j = 0; j < tmpNeighbors.length; j++) {
				neighbor = tmpNeighbors[j];
				if (neighbor.visitId !== idVisit) {
					neighbor.visitId = idVisit;
					dist = distPlane(point, neighbor);
					if (dist > 0) {
						neighbor.visible = true;
						visible.push(neighbor);
						
					}
				}
			}
		}
		
		return visible;
	};
	
	/**
	 * return vertices for intersection of two simplicial facets
	 * 
	 * @param {Facet} facetA - first facet
	 * @param {Facet} facetB - second facet
	 * @return {Array|Float32Array} the vertices of interserction
	 * 
	 */ 
	
	var facetIntersect = quickhull._utils.facetIntersect = function(facetA, facetB) {
		var intersection = new Array();
		var verticesA = facetA.vertices;
		var verticesB = facetB.vertices;
		var tmpVertex;
		
		for ( var i = 0; i < verticesA.length; i++) {
			tmpVertex = verticesA[i];
			
			if(isIn(verticesB, tmpVertex)){
				intersection.push(tmpVertex);
			}
		}
		
		return intersection;
	};
	
	/**
	 * 
	 * 
	 */
	
	var orientNewFacet = quickhull._utils.orientFacet = function(visible, newFacet, intersection){
		var interiorPoint = getCenter(visible.vertices);
		
		setFacetPlane(newFacet);
		
		var dist = distPlane(interiorPoint, newFacet);
		
		if (dist > 0 ) {
			newFacet.offset = -newFacet.offset;
			newFacet.topoOriented = false;
			for ( var k = 0; k < newFacet.normal.length; k++) {
				newFacet.normal[k] = -newFacet.normal[k];
			}
		}
	};
	
	/**
	 * make new facets for simplicial visible facet and apex
	 * 
	 * @param {Array} point - the point
	 * @param {Facet} visible - facet visible form point
	 * @return {Array|Facet} the new facets
	 * 
	 */ 
	
	var makeNewSimplicial = quickhull._utils.makeNewSimplicial = function(apex, visible) {
		var tmpNeighbors = visible.neighbors;
		var neighbor, newFacet;
		var apexArr = new Array(), tmpVertices, vertices, newFacets = new Array();
		apexArr.push(apex);
				
		for ( var i = 0; i < tmpNeighbors.length; i++) {
			neighbor = tmpNeighbors[i];
			if (!(neighbor.visible)) {
				tmpVertices = facetIntersect(neighbor, visible);
				
				vertices = apexArr.concat(tmpVertices);

				newFacet = new Facet(); 
				newFacet.vertices = vertices;
				newFacet.neighbors.push(neighbor);

				orientNewFacet(visible, newFacet, visible.vertices);

				neighbor.neighbors.push(newFacet);

				newFacets.push(newFacet);
			}
		}
		
		return newFacets;
	};
	
	/**
	 * make new facets from point
	 * 
	 * @param {Array} point - the point
	 * @param {Array|Facet} visible - facets visible form point
	 * @return {Array|Facet} the new facets
	 * 
	 */ 
	
	var makeNewFacets = quickhull._utils.makeNewFacets = function(point, visible) {
		
		var newFacets = new Array(), tmpFacets;
		
		for ( var i = 0; i < visible.length; i++) {
			tmpFacets = makeNewSimplicial(point, visible[i]);
			
			newFacets = newFacets.concat(tmpFacets);
		}
		
		return newFacets;
	};
	
	/**
	 * match newfacets to their newfacet neighbors
	 * 
	 * @param {Array|Facet} newFacets - the new facets
	 * 
	 */ 
	
	var matchNewFacets = quickhull._utils.matchNewFacets = function(newFacets) {
		var intersect;
		var numVert = newFacets[0].vertices.length;
		
		for ( var i = 0; i < newFacets.length; i++) {
			for ( var j = 0; j < newFacets.length; j++) {
				if (j !== i) {
					intersect = facetIntersect(newFacets[i], newFacets[j]);
					
					if (intersect.length === numVert -1) {
						newFacets[i].neighbors.push(newFacets[j]);
					}
				}
			}
		}
	};
	
	/**
	 * update vertices' set
	 * 
	 * @param {Array|Float32Array} currentVertices - the current set of vertices
	 * @param {Array|Facet} visible - the facets visible from new vertex
	 * @param {Array|Facet} newFacets - the new facets
	 * 
	 */ 
	
	var updateVertices = quickhull._utils.updateVertices = function(currentVertices, visible, newFacets) {
		var visibleVertices = new Array();
		var newFacetVertices = new Array();
		
		for ( var i = 0; i < visible.length; i++) {
			for ( var j = 0; j < visible[i].vertices.length; j++) {
				addNoDup(visibleVertices, visible[i].vertices[j]);
			}
		}
		
		for ( var i = 0; i < newFacets.length; i++) {
			for ( var j = 0; j < newFacets[i].vertices.length; j++) {
				addNoDup(newFacetVertices, newFacets[i].vertices[j]);
			}
		}
		
		for ( var i = 0; i < visibleVertices.length; i++) {
			if (!isIn(newFacetVertices, visibleVertices[i])) {
				removePoint(currentVertices, visibleVertices[i]);
			}
		}
	};
	
	/**
	 * partitions points in visible facets to new facets
	 * 
	 * @param {Array|Facet} visible - the facets visible from new vertex
	 * @param {Array|Facet} newFacets - the new facets
	 * @param {Array} newVertex - the new vertex added to the hull
	 * 
	 */ 
	
	var partitionVisible = quickhull._utils.partitionVisible = function(visible, newFacets, newVertex) {
		// insieme di tutti i punti esterni, senza ripetizioni, delle facce visibili dal nuovo vertice inserito
		var outsidePoints = new Array();
		var tmpOutside;
		
		for ( var i = 0; i < visible.length; i++) {
			if (visible[i].furthestPoint !== null && (comparePoints(newVertex, visible[i].furthestPoint) !== 0)) {
				addNoDup(outsidePoints, visible[i].furthestPoint);
			}
			
			tmpOutside = visible[i].outsideSet;
			
			for ( var j = 0; j < tmpOutside.length; j++) {
				if (comparePoints(newVertex, tmpOutside[j]) !== 0) {
					addNoDup(outsidePoints, tmpOutside[j]);
				}
			}
		}
		
		for ( var i = 0; i < newFacets.length; i++) {
			for ( var j = 0; j < outsidePoints.length; j++) {
				addOutside(outsidePoints[j], newFacets[i]);
			}
		}
	};
	
	/**
	 * remove visible facets from current convex hull
	 * 
	 * @param {Array|Facet} currentHull - the facet to add point
	 * @param {Array|Facet} visible - set of visible facet
	 */ 
	
	var deleteVisible = quickhull._utils.deleteVisible = function(currentHull, visible) {
		var tmpNeighbors;
		
		for ( var i = 0; i < currentHull.length; i++) {
			if (currentHull[i].visible) {// se è visibile viene eliminta dal guscio
				currentHull.splice(i,1);
				i = i -1;
			} else {// altrimenti vengono eliminati tutti i vicini visibili
				tmpNeighbors = currentHull[i].neighbors;
				for ( var j = 0; j < tmpNeighbors.length; j++) {
					if (tmpNeighbors[j].visible) {// se è visibile viene eliminta dall'insieme dei vicini
						tmpNeighbors.splice(j,1);
						j = j -1;
					}
				}
			}
		}
		
		
	};
	
	/**
	 * add point (usually furthest point) above facet to hull
	 * 
	 * @param {Facet} facet - the facet to add point
	 * @param {Array|Facet} currentHull - the current hull
	 * @param {Array|Float32Array} currentVertices - the current vertices of the current hull
	 * 
	 */ 
	
	var addPoint = quickhull._utils.addPoint = function(facet, currentHull, currentVertices) {
		var vertex = facet.furthestPoint;
		
		var horizon = findHorizon(vertex, facet);
		
		if (displayGlob === 2) {
			for ( var i = 1; i < horizon.length; i++) {
				//simplPlasm = SIMPLICIAL_COMPLEX(horizon[i].vertices)([[0,1,2]]);
				simplPlasm = SIMPLICIAL_COMPLEX(horizon[i].vertices)(cellsPlasm(horizon[i].vertices,[horizon[i]]));
				steps.push(STRUCT([steps[steps.length - 1], COLOR([0,0,0])(SKELETON(1)(simplPlasm)), COLOR([1,0,0,0.8])(simplPlasm)]));
			}
		}
		
		var newFacets = makeNewFacets(vertex, horizon);
		
		if (displayGlob === 2) {
			for ( var i = 0; i < newFacets.length; i++) {
				//simplPlasm = SIMPLICIAL_COMPLEX(newFacets[i].vertices)([[0,1,2]]);
				simplPlasm = SIMPLICIAL_COMPLEX(newFacets[i].vertices)(cellsPlasm(newFacets[i].vertices,[newFacets[i]]));
				steps.push(STRUCT([steps[steps.length - 1], COLOR([0,0,0])(SKELETON(1)(simplPlasm)), COLOR([0,0,1,0.8])(simplPlasm)]));
			}
		}
		
		if (displayGlob === 1 || displayGlob === 3 ) {
			// disegno delle nuove facce su canvas
			for ( var i = 0; i < newFacets.length; i++) { 
				plotBaseLine(newFacets[i].vertices,'rgb(0,0,180)');
			}
		}
		
		matchNewFacets(newFacets);
		
		currentVertices.push(vertex);
		
		updateVertices(currentVertices, horizon, newFacets);
		
		partitionVisible(horizon, newFacets, vertex);
		
		
		deleteVisible(currentHull);
		
		//aggiunge le nuove facce al guscio
		for ( var i = 0; i < newFacets.length; i++) {
			currentHull.push(newFacets[i]);
		}
	};
	
	/**
	 * calculate the cell's facet for plasm's SIMPLICIAL_COMPLEX from global vertices
	 * 
	 * @param {Array|Float32Array} vertices - the current vertices
	 * @param {Array|Facet} visible - set of visible facets
	 * @return {Array|Float32Array} the vertices global index of facet
	 */ 
	
	var cellsPlasm = quickhull._utils.cellPlasm = function(vertices, facets) {
		var tmpVertices;
		var cells = new Array();
		for ( var i = 0; i < facets.length; i++) {
			tmpVertices = facets[i].vertices;
			cells.push([]);
			for ( var j = 0; j < tmpVertices.length; j++) {
				for ( var k = 0; k < vertices.length; k++) {
					if ((comparePoints(vertices[k], tmpVertices[j]) === 0)) {
						cells[i].push(k);
					}
				}
			}
		}

		return cells;
	};
	
	/**
	 * Algorithm for convex hull computing
	 * 
	 * @param {Array|Float32Array} points - set of points with the same dimension
	 * @param {Number} display - 0|1|2|3 number for draw convex hull: 0 nothing, 1 canvas, 2 plasm step by step, 3 canvas & plasm
	 * @param {Number} log - 0|1 number for log print
	 * @return {Array|Facet} convexHull
	 * 
	 */ 

	var quickhull = quickhull.quickhull = function(points, display, log) {
		//TODO verifica punti tutti della stessa dimensione
		var numPoints = points.length;
		var dim = points[0].length;

		var maximum, minimum, maxsMins = new Array();
		var point;
		var vertices = new Array();
		var convexHull = new Array();
		var maxCoord = -Infinity, minCoord = Infinity, maxDet, det;
		var maxX, minX, maxPoint;
		var nextFacet;
		var simplPlasm;
		
		if (display !== undefined) {
			displayGlob = display;
		}
		
		if (log !== undefined) {
			logGlob = log;
		}
		
		if (displayGlob === 2) {
			steps = new Array();
		}
		
		// costruzione insiemi max min per ogni dimensione
		for ( var i = 0; i < dim; i++) {
			maximum = minimum = points[0];
			for ( var j = 0; j < numPoints; j++) {
				point = points[j];
				if (maximum[i] < point[i]) {
					maximum = point;
				} else if (minimum[i] > point[i]) {
					minimum = point;
				}
			}
			maxsMins.push(maximum);
			maxsMins.push(minimum);
		}

		// vertici iniziali del simplesso/guscio iniziale composto da dim + 1 punti
		if (maxsMins.length >= 2) {
			for ( var j = 0; j < maxsMins.length; j++) {
				point = maxsMins[j];
				if (maxCoord < point[0]) {
					maxCoord = point[0];
					maxX = point;
				} else if (minCoord > point[0]) {
					minCoord = point[0];
					minX = point;
				}
			}
		}else{
			for ( var j = 0; j < numPoints; j++) {
				point = points[j];
				if (maxCoord < point[0]) {
					maxCoord = point[0];
					maxX = point;
				} else if (minCoord > point[0]) {
					minCoord = point[0];
					minX = point;
				}
			}
		}

		addNoDup(vertices, maxX);
		addNoDup(vertices, minX);
		
		for ( var k = vertices.length; k < dim + 1; k++) {
			maxPoint = null;
			maxDet = -Infinity;
			for ( var j = 0; j < maxsMins.length; j++) {
				point = maxsMins[j];
				if (!isIn(vertices, point)) {
					det = detSimplex(point, vertices, k); 
					det = Math.abs(det); 
					if (det > maxDet) {
						maxDet = det;
						maxPoint = point;
						//TODO gestire determinanti molto vicini a zero (nearzero)
					}
				}
			}

			if (maxPoint === null) {
				for ( var j = 0; j < numPoints; j++) {
					point = points[j];
					if (!isIn(vertices, point)) {
						det = detSimplex(point, vertices, k); 
						det = Math.abs(det); 
						if (det > maxDet) {
							maxDet = det;
							maxPoint = point;
						}
					}
				}
			}

			//TODO gestire il caso in cui maxPoint non viene trovato
			if (maxPoint !== null) {
				addNoDup(vertices, maxPoint);
			}
		}
		
		if (logGlob === 1) {
			console.log("Simplesso/vertici iniziale/i: ");
			for ( var i = 0; i < vertices.length; i++) {
				console.log("[" + vertices[i] + "]");
			}
		}
		
		convexHull = initialHull(vertices);
		
		if (displayGlob === 1 || displayGlob === 3 ) {
			// disegno su canvas del guscio iniziale
			for ( var i = 0; i < convexHull.length; i++) { 
				plotBaseLine(convexHull[i].vertices,'rgb(180,180,180)');
			}
		}
		
		partitionAll(convexHull, vertices, points);
		
		nextFacet = furthestNext(convexHull);
		
		while (nextFacet !== null) {
			
			if (displayGlob === 2) {
				simplPlasm = SIMPLICIAL_COMPLEX(vertices)(cellsPlasm(vertices,convexHull));
				steps.push(STRUCT([COLOR([1,1,1,0.5])(simplPlasm), COLOR([0,0,0])(SKELETON(1)(simplPlasm))]));
			}
			
			if (displayGlob === 2) {
				//simplPlasm = SIMPLICIAL_COMPLEX(nextFacet.vertices)([[0,1,2]]);
				simplPlasm = SIMPLICIAL_COMPLEX(nextFacet.vertices)(cellsPlasm(nextFacet.vertices,[nextFacet]));
				steps.push(STRUCT([steps[steps.length - 1], COLOR([0,0,0])(SKELETON(1)(simplPlasm)), COLOR([1,0,0,0.8])(simplPlasm)]));
			}
			
			if (displayGlob === 1 || displayGlob === 3 ) {
				// disegno del guscio corrente su canvas e dalla faccia con il punto più lontano
				qhPlotPoints(points);

				for ( var i = 0; i < convexHull.length; i++) {
					plotBaseLine(convexHull[i].vertices,'rgb(180,180,180)');
				}
				
				plotBaseLine(nextFacet.vertices,'rgb(255,0,0)');
			}
			
			addPoint(nextFacet, convexHull, vertices);
			
			nextFacet = furthestNext(convexHull);
			
			
		}
		
		if (logGlob === 1) {
			console.log("Vertici: ");
			console.log(vertices);
		}
		
		if (displayGlob === 2) {
			console.log(cellsPlasm(vertices,convexHull));
			simplPlasm = SIMPLICIAL_COMPLEX(vertices)(cellsPlasm(vertices,convexHull));
			steps.push(STRUCT([COLOR([0,1,0,0.5])(simplPlasm), COLOR([0,0,0])(SKELETON(1)(simplPlasm))]));
		}
		
		return {vertices:vertices,facets:convexHull,steps:steps};
	};
	
	

}(this));