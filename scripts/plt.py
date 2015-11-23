import cairo

w = 800
h = 600
p = { 'n': 50, 's': 40, 'e': 20, 'w': 50 }

ims = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
cr = cairo.Context(ims)

# outline
cr.set_source_rgb(1, 1, 1)
cr.rectangle(0, 0, w, h)
cr.fill()
cr.set_source_rgb(0, 0, 0)
cr.rectangle(p['n'], p['w'], w - p['e'] - p['w'], h - p['s'] - p['n'])
cr.stroke()

# title
cr.set_source_rgb(0, 0, 0)
cr.select_font_face("Open Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
cr.set_font_size(20)
cr.move_to(p['w'], p['n'] / 1.5)
cr.show_text("Convergence of k_on for Sialic Acid to avian N1")

ims.write_to_png("image.png")

