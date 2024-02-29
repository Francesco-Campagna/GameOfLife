#ifndef JALLEGRO_GRAPHICS_H
#define JALLEGRO_GRAPHICS_H

#include <allegro5/allegro_primitives.h>
#include "Color.h"

namespace jallegro {
    class Graphics {

    public:
        void fillRect(int x, int y, int w, int h) {
            al_draw_filled_rectangle(x, y, x + w, y + h, color.toAllegro());
        }

        void fillCircle(int x, int y, int r) {
            al_draw_filled_circle(x, y, r, color.toAllegro());
        }

        void setColor(Color c) {
            color = c;
        }

    private:
        Color color;
    };
}

#endif