/*
 *Copyright © Eliott Veyrier
 *Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *The Software is provided “as is”, without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders X be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the Software. »
*/


shader_type particles;
uniform vec3 gravity = vec3(0,-9.8,0);
uniform vec3 initial_velocity = vec3(5,10,0);
uniform float initial_velocity_scale = 1.;
uniform float initial_velocity_randomness = 1.;
uniform float drag : hint_range(0.5,5.5,0.1)= 3.;
uniform float scatter : hint_range(0.,2.,0.1)= 1.;
//this gradient gets sampled randomly
uniform sampler2D random_colors : hint_black;
//rng from the ParticleMaterial
float rand_from_seed(in uint seed) {
  int k;
  int s = int(seed);
  if (s == 0)
    s = 305420679;
  k = s / 127773;
  s = 16807 * (s - k * 127773) - 2836 * k;
  if (s < 0)
    s += 2147483647;
  seed = uint(s);
  return float(seed % uint(65536)) / 65535.0;
}
uint hash(uint x) {
  x = ((x >> uint(16)) ^ x) * uint(73244475);
  x = ((x >> uint(16)) ^ x) * uint(73244475);
  x = (x >> uint(16)) ^ x;
  return x;
}
/* simplex noise under The MIT License
 * Copyright © 2013 Nikita Miropolskiy */
/* discontinuous pseudorandom uniformly distributed in [-0.5, +0.5]^3 */
vec3 random3(vec3 c) {
	float j = 4096.0*sin(dot(c,vec3(17.0, 59.4, 15.0)));
	vec3 r;
	r.z = fract(512.0*j);
	j *= .125;
	r.x = fract(512.0*j);
	j *= .125;
	r.y = fract(512.0*j);
	return r-0.5;
}
/* skew constants for 3d simplex functions */
const float F3 =  0.3333333;
const float G3 =  0.1666667;
/* 3d simplex noise */
float simplex3d(vec3 p) {
	 /* 1. find current tetrahedron T and it's four vertices */
	 /* s, s+i1, s+i2, s+1.0 - absolute skewed (integer) coordinates of T vertices */
	 /* x, x1, x2, x3 - unskewed coordinates of p relative to each of T vertices*/
	 
	 /* calculate s and x */
	 vec3 s = floor(p + dot(p, vec3(F3)));
	 vec3 x = p - s + dot(s, vec3(G3));
	 
	 /* calculate i1 and i2 */
	 vec3 e = step(vec3(0.0), x - x.yzx);
	 vec3 i1 = e*(1.0 - e.zxy);
	 vec3 i2 = 1.0 - e.zxy*(1.0 - e);
	 	
	 /* x1, x2, x3 */
	 vec3 x1 = x - i1 + G3;
	 vec3 x2 = x - i2 + 2.0*G3;
	 vec3 x3 = x - 1.0 + 3.0*G3;
	 
	 /* 2. find four surflets and store them in d */
	 vec4 w, d;
	 
	 /* calculate surflet weights */
	 w.x = dot(x, x);
	 w.y = dot(x1, x1);
	 w.z = dot(x2, x2);
	 w.w = dot(x3, x3);
	 
	 /* w fades from 0.6 at the center of the surflet to 0.0 at the margin */
	 w = max(0.6 - w, 0.0);
	 
	 /* calculate surflet components */
	 d.x = dot(random3(s), x);
	 d.y = dot(random3(s + i1), x1);
	 d.z = dot(random3(s + i2), x2);
	 d.w = dot(random3(s + 1.0), x3);
	 
	 /* multiply d by w^4 */
	 w *= w;
	 w *= w;
	 d *= w;
	 
	 /* 3. return the sum of the four surflets */
	 return dot(d, vec4(52.0));
}
float simplex3d_fractal(vec3 m) {
    return   0.5333333*simplex3d(m)
			+0.2666667*simplex3d(2.0*m)
			+0.1333333*simplex3d(4.0*m);
}

/**/
//Rotation around an axis
mat4 rotationAxisAngle( vec3 v, float angle ){
    float s = sin( angle );
    float c = cos( angle );
    float ic = 1.0 - c;
	return mat4(
		vec4(v.x*v.x*ic + c,v.x*v.y*ic + s*v.z,v.x*v.z*ic - s*v.y,0.),
		vec4(v.y*v.x*ic - s*v.z,v.y*v.y*ic + c,v.y*v.z*ic + s*v.x,0.0),
		vec4(v.z*v.x*ic + s*v.y,v.z*v.y*ic - s*v.x,v.z*v.z*ic + c,0.0),
		vec4(0.0,0.0,0.0,1.0)
		);
}
//    return mat4(v.x*v.x*ic + c/**/, v.y*v.x*ic - s*v.z, v.z*v.x*ic + s*v.y, 0.0,
//                v.x*v.y*ic + s*v.z, v.y*v.y*ic + c,     v.z*v.y*ic - s*v.x, 0.0,
//                v.x*v.z*ic - s*v.y, v.y*v.z*ic + s*v.x, v.z*v.z*ic + c,     0.0,
//		          0.0,                0.0,                0.0,                1.0 );

//box muhler transform to get gaussian distributed 2D point
vec2 rand_box_muhler(uint h1,uint h2){
	vec2 U = vec2(rand_from_seed(h1),rand_from_seed(h2));
	return vec2(
		sqrt(-2.*log(U.x))*cos(U.y*6.28318530718),
		sqrt(-2.*log(U.y))*sin(U.x*6.28318530718)
		);
}
void vertex(){
	if (RESTART){
		//seeds for the rng
		uint alt_seed1 = hash(NUMBER + uint(1) + RANDOM_SEED);
		uint alt_seed2 = hash(NUMBER + uint(27) + RANDOM_SEED);
		uint alt_seed3 = hash(NUMBER + uint(43) + RANDOM_SEED);
		uint alt_seed4 = hash(NUMBER + uint(111) + RANDOM_SEED);
		VELOCITY = initial_velocity_scale *( 
			initial_velocity
			 + initial_velocity_randomness*rand_box_muhler(alt_seed2,alt_seed3).xyx
			);
		COLOR = texture(random_colors,vec2(rand_from_seed(alt_seed1),0));
	} else {
	vec3 normal = normalize(TRANSFORM[2].xyz);
	// smooth random rotation that cheaply emulates turbulence rotating the confetti
	TRANSFORM *= rotationAxisAngle(normalize(VELOCITY), 10.*DELTA *simplex3d_fractal(TRANSFORM[3].xyz));
	TRANSFORM *= rotationAxisAngle(normal, 10.*DELTA *simplex3d_fractal(float(NUMBER) + TRANSFORM[3].xyz));
	CUSTOM.y += DELTA;
	// gravity and drag, drag is proportional to how much confettis are facing the wind hence the dot product
	VELOCITY += DELTA * gravity - DELTA * drag * abs(dot(VELOCITY,normal))*normalize(VELOCITY) ;
	VELOCITY += scatter * DELTA * (TRANSFORM * vec4(VELOCITY,0.)).xyz
	}
	
	
}
