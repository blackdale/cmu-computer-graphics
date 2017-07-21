#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i], Matrix3x3::identity());
  }

  // set top level transformation
  transformation = canvas_to_screen;

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y++;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y++;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y--;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y--;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::setup_supersample()
{
	delete this->supersample_render_target;
	this->supersample_render_target = 0;
	if (sample_rate > 1) {
		this->supersample_render_target = new unsigned char[4 * sample_rate * sample_rate * target_w * target_h];
		memset(supersample_render_target, 255, 4 * target_w * target_h * sample_rate * sample_rate);
	}
}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  setup_supersample();
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  setup_supersample();
}

//Matrix3x3 groupTransform = Matrix3x3::identity();

void SoftwareRendererImp::draw_element( SVGElement* element, Matrix3x3 groupTransform ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack


	if (element->type != GROUP) {
		transformation = canvas_to_screen * groupTransform * element->transform;
	}


  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element), groupTransform * element->transform);
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );
}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group, Matrix3x3 groupTransform ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i], groupTransform);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  unsigned char* rt = sample_rate > 1 ? supersample_render_target : render_target;
  int w = target_w * sample_rate;
  int h = target_h * sample_rate;

  rasterize_pointi(rt, w, h, sx, sy, color);
}

void SoftwareRendererImp::rasterize_pointi(unsigned char* rt, int w, int h, int sx, int sy, Color color) {
	// check bounds
	if (sx < 0 || sx >= w) return;
	if (sy < 0 || sy >= h) return;

	//unsigned char* rt = sample_rate > 1 ? supersample_render_target : render_target;
	//unsigned char* rt = render_target;

	// fill sample - NOT doing alpha blending!
	rt[4 * (sx + sy * w)] = (uint8_t)(color.r * 255);
	rt[4 * (sx + sy * w) + 1] = (uint8_t)(color.g * 255);
	rt[4 * (sx + sy * w) + 2] = (uint8_t)(color.b * 255);
	rt[4 * (sx + sy * w) + 3] = (uint8_t)(color.a * 255);
}

void SoftwareRendererImp::rasterize_line( float x0f, float y0f,
                                          float x1f, float y1f,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization

	int x0 = floor(x0f * sample_rate);
	int x1 = floor(x1f * sample_rate);
	int y0 = floor(y0f * sample_rate);
	int y1 = floor(y1f * sample_rate);

	bool steep = abs(x1 - x0) > abs(y1 - y0);
	if (!steep) {
		swap(x0, y0);
		swap(x1, y1);
	}
	if (x0 > x1) {
		swap(x0, x1);
		swap(y0, y1);
	}

	int dx = x1 - x0,
		dy = y1 - y0,
		y = y0,
		eps = 0;

	int diff = dy > 0 ? 1 : -1;

	unsigned char* rt = sample_rate > 1 ? supersample_render_target : render_target;
	int w = target_w * sample_rate;
	int h = target_h * sample_rate;
	
	for (int x = x0; x <= x1; x++) {
		if (steep)
			rasterize_pointi(rt, w, h, x, y, color);
		else
			rasterize_pointi(rt, w, h, y, x, color);

		eps += dy;
		if (diff * (eps << 1) >= dx * diff) {
			y += diff;  eps -= dx * diff;
		}
	}
}

enum Status {
	FullyAccept,
	Accept,
	Reject
};

bool is_trivial_reject_point(float a, float b, float c, float x, float y) {
	return a * x + b * y + c < 0;
}

Status is_trivial_accept_tile(int tile_x, int tile_y, int tile_w, int tile_h, float a, float b, float c) {
	bool x0y0 = !is_trivial_reject_point(a, b, c, tile_x, tile_y);
	bool x1y0 = !is_trivial_reject_point(a, b, c, tile_x + tile_w, tile_y);
	bool x0y1 = !is_trivial_reject_point(a, b, c, tile_x, tile_y + tile_h);
	bool x1y1 = !is_trivial_reject_point(a, b, c, tile_x + tile_w, tile_y + tile_h);

	if (x0y0 && x1y0 && x0y1 && x1y1)
		return FullyAccept;
	return Accept;
}

bool is_trivial_reject_tile(int tile_x, int tile_y, int tile_w, int tile_h, float a, float b, float c) {
	return is_trivial_reject_point(a, b, c, tile_x, tile_y) &&
		is_trivial_reject_point(a, b, c, tile_x + tile_w, tile_y) &&
		is_trivial_reject_point(a, b, c, tile_x, tile_y + tile_h) &&
		is_trivial_reject_point(a, b, c, tile_x + tile_w, tile_y + tile_h);
}

void SoftwareRendererImp::rasterize_triangle( float x0i, float y0i,
                                              float x1i, float y1i,
                                              float x2i, float y2i,
                                              Color color ) {
  // Task 3: 
  // Implement triangle rasterization
  
	/*rasterize_line(x0, y0, x1, y1, color);
	rasterize_line(x1, y1, x2, y2, color);
	rasterize_line(x2, y2, x0, y0, color);*/

	int tile_w = 16;
	int tile_h = 16;

	unsigned char* rt = (sample_rate > 1) ? supersample_render_target : render_target;
	int w = target_w * sample_rate;
	int h = target_h * sample_rate;

	int tile_count_w = w / tile_w + 1;
	int tile_count_h = h / tile_h + 1;

	std::vector<std::vector<bool>> tiles(tile_count_w);
	for (int i = 0; i < tile_count_h; ++i) {
		tiles.at(i).reserve(tile_count_h);
	}

	float x0 = x0i * sample_rate;
	float x1 = x1i * sample_rate;
	float x2 = x2i * sample_rate;

	float y0 = y0i * sample_rate;
	float y1 = y1i * sample_rate;
	float y2 = y2i * sample_rate;

	float a0 = y0 - y1;
	float b0 = x1 - x0;
	float c0 = y1 * x0 - y0 * x1;

	float a1 = y1 - y2;
	float b1 = x2 - x1;
	float c1 = y2 * x1 - y1 * x2;

	float a2 = y2 - y0;
	float b2 = x0 - x2;
	float c2 = y0 * x2 - y2 * x0;

	float det = b0 * a2 - b2 * a0;

	if (det < 0) {
		a0 = -a0; a1 = -a1; a2 = -a2;
		b0 = -b0; b1 = -b1; b2 = -b2;
		c0 = -c0; c1 = -c1; c2 = -c2;
	}
	
	for (int i = 0; i < tile_count_w; ++i) {
		for (int j = 0; j < tile_count_h; ++j) {
			bool trivial_reject = is_trivial_reject_tile(i * tile_w, j * tile_h, tile_w, tile_h, a0, b0, c0) ||
				is_trivial_reject_tile(i * tile_w + 0.5f, j * tile_h + 0.5f, tile_w, tile_h, a1, b1, c1) ||
				is_trivial_reject_tile(i * tile_w + 0.5f, j * tile_h + 0.5f, tile_w, tile_h, a2, b2, c2);

			if (!trivial_reject) {
				Status s0 = is_trivial_accept_tile(i * tile_w + 0.5f, j * tile_h + 0.5f, tile_w, tile_h, a0, b0, c0);
				Status s1 = is_trivial_accept_tile(i * tile_w + 0.5f, j * tile_h + 0.5f, tile_w, tile_h, a1, b1, c1);
				Status s2 = is_trivial_accept_tile(i * tile_w + 0.5f, j * tile_h + 0.5f, tile_w, tile_h, a2, b2, c2);

				if (s0 == FullyAccept && s1 == FullyAccept && s2 == FullyAccept) {
					for (int x = 0; x < tile_w; ++x) {
						for (int y = 0; y < tile_h; ++y) {
							rasterize_pointi(rt, w, h, i * tile_w + x, j * tile_h + y, color);
						}
					}
				}
				else if (s0 == Accept || s1 == Accept || s2 == Accept) {
					for (int x = 0; x < tile_w; ++x) {
						for (int y = 0; y < tile_h; ++y) {
							if (!is_trivial_reject_point(a0, b0, c0, i * tile_w + x + 0.5f, j * tile_h + y + 0.5f) &&
								!is_trivial_reject_point(a1, b1, c1, i * tile_w + x + 0.5f, j * tile_h + y + 0.5f) &&
								!is_trivial_reject_point(a2, b2, c2, i * tile_w + x + 0.5f, j * tile_h + y + 0.5f))
							rasterize_pointi(rt, w, h, i * tile_w + x, j * tile_h + y, color);
						}
					}
				}
			}
		}
	}
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
	
	if (sample_rate <= 1)
		return;

	int w = target_w * sample_rate;
	int h = target_h * sample_rate;
	for (int i = 0; i < w; i += sample_rate) {
		for (int j = 0; j < h; j += sample_rate) {
			float r = 0, g = 0, b = 0, a = 0;
			for (int x = 0; x < sample_rate; ++x) {
				for (int y = 0; y < sample_rate; ++y) {
					r += (float)supersample_render_target[4 * (i + x + (j + y) * w)] / 255.f;
					g += (float)supersample_render_target[4 * (i + x + (j + y) * w) + 1] / 255.f;
					b += (float)supersample_render_target[4 * (i + x + (j + y) * w) + 2] / 255.f;
					a += (float)supersample_render_target[4 * (i + x + (j + y) * w) + 3] / 255.f;
				}
			}
			r /= sample_rate * sample_rate;
			g /= sample_rate * sample_rate;
			b /= sample_rate * sample_rate;
			a /= sample_rate * sample_rate;

			rasterize_pointi(render_target, target_w, target_h, i / sample_rate, j / sample_rate, Color(r, g, b, a));
		}
	}
	
  return;
}

} // namespace CMU462
