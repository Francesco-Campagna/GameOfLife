#ifndef JALLEGRO_COLOR_H
#define JALLEGRO_COLOR_H

namespace jallegro {
    class Color {

    public:
        Color() {
            color = al_map_rgb(0, 0, 0);
        }

        Color(unsigned int r, unsigned int g, unsigned int b) {
            color = al_map_rgb(0, 0, 0);
            if (r > 255 || g > 255 || b > 255)
                return;
            color = al_map_rgb(r, g, b);
        }

        static Color blue() {
            return Color(0, 0, 255);
        }

        static Color red() {
            return Color(255, 0, 0);
        }

        static Color green() {
            return Color(0, 255, 0);
        }

        static Color black() {
            return Color(0, 0, 0);
        }

        static Color white() {
            return Color(255, 255, 255);
        }

        static Color yellow() {
            return Color(255, 255, 0);
        }

        ALLEGRO_COLOR toAllegro() {
            return color;
        }

    private:
        ALLEGRO_COLOR color;

    };
}
#endif