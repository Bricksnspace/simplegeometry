/*
	Copyright 2014 Mario Pascucci <mpascucci@gmail.com>
	This file is part of JSimple3DGeom

	JSimple3DGeom is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JSimple3DGeom is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JSimple3DGeom.  If not, see <http://www.gnu.org/licenses/>.

*/


package bricksnspace.j3dgeom;


/**
 * A simple 3D point
 * 
 * @author Mario Pascucci
 */
public class Point3D {
	
	public static final Point3D ORIGIN = new Point3D(0,0,0); 

	public float x,y,z;
	
	
	public Point3D(float px, float py, float pz) {
		
		x = px;
		y = py;
		z = pz;
	}
	
	
	
	public Point3D(Point3D p) {
		
		x = p.x;
		y = p.y;
		z = p.z;		
	}
	
	
	public Point3D(float[] p) {
		
		x = p[0];
		y = p[1];
		z = p[2];
	}
	
	
	
	public Point3D(float[] p, int offset) {
		
		x = p[offset];
		y = p[offset+1];
		z = p[offset+2];
	}
	
	
	@Override
	public String toString() {
//		return String.format("Pt3D [x=%s, y=%s, z=%s]", x, y, z);
		return String.format("Pt3D [x=%.3f, y=%.3f, z=%.3f]", x, y, z);

	}



	public Point3D transform(Matrix3D m) {
		
		return new Point3D(m.transformPoint(x, y, z));
	}
	
	
	
	public float getDistSq(Point3D p) {
		
		float dx = x-p.x;
		float dy = y-p.y;
		float dz = z-p.z;
		return (dx*dx+dy*dy+dz*dz);
	}
	
	
	
	public Point3D translate(Point3D offset) {
		
		return new Point3D(x+offset.x,y+offset.y,z+offset.z);
	}
	
	
	
	public boolean coincident(Point3D p) {
		
		return x == p.x && y == p.y && z == p.z;
	}
	
	
	/**
	 * Return a vector oriented from this point to point p
	 * @param p head of vector
	 * @return a vector
	 */
	public Point3D vector(Point3D p) {
		
		return new Point3D(p.x-x,p.y-y,p.z-z);
	}
	
	
	/**
	 * returns modulo of vector as a point
	 */
	public float modulo() {
		return (float) Math.sqrt(x*x+y*y+z*z);
	}
	
	
	public Point3D normalize() {
		
		float mod = modulo();
		if (mod == 0) {
			return new Point3D(0,0,0);
		}
		return new Point3D(x/mod,y/mod,z/mod);
	}



	public Point3D scale(float factor) {

		return new Point3D(x*factor, y*factor, z*factor);
	}
	
}
