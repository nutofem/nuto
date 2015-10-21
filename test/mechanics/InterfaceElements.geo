
meshSize    = 2;
height      = 10;
width       = 10;

///////////////////////////////
//   POINTS                  //
///////////////////////////////

// RECTANGLE
Point(1)  = {0    	  ,   0    		    ,   0,  meshSize}; 
Point(2)  = {width    ,   0    		    ,   0,  meshSize}; 
Point(3)  = {width    ,   height/2. 	,   0,  meshSize}; 
Point(4)  = {width    ,   height    	,   0,  meshSize}; 
Point(5)  = {0  	  ,   height     	,   0,  meshSize}; 
Point(6)  = {0   	  ,	  height/2.    	,   0,  meshSize}; 

// LINE
Point(7)  = {width/2  ,   height/2.     ,   0,  meshSize};

///////////////////////////////
//   LINES                   //
///////////////////////////////

// RECTANGLE
Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};
Line(5)  = {5,6};
Line(6)  = {6,1};

// LINE
Line(7)  = {7,3};

///////////////////////////////
//   LOOPS                   //
///////////////////////////////

Line Loop(1) = {1:6};

///////////////////////////////
//   PLANES                  //
///////////////////////////////

Plane Surface(1) = {1};
Line {7} In Surface {1};

///////////////////////////////
//   PHYSICAL GROUPS         //
///////////////////////////////

Physical Line(777)    = {7};
Physical Surface(999) = {1};



