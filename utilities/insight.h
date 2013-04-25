/*
 * 07/11
 *
 * Insight3 file format constants.
 *
 *
 * Hazen
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc tracker.c -o tracker
 *
 */

// Header
#define VERSION 0
#define FRAMES 4
#define STATUS 8
#define MOLECULES 12
#define DATA 16

// Localizations
#define OBJECT_DATA_SIZE 18
#define DATUM_SIZE 4

#define XO 0        // original x coordinate
#define YO 1        // original y coordinate
#define X 2         // drift corrected x coordinate
#define Y 3         // drift corrected y coordinate

#define HEIGHT 4    // fit peak height
#define AREA 5      // fit peak area
#define WIDTH 6     // fit peak width
#define ANGLE 7     // fit peak angle (unused, so repurposed as a visited flag)
#define VISITED 7   // visited flag

#define ASPECT 8    // fit peak aspect ratio
#define BACKG 9     // fit peak background
#define SUM 10      // sum - background for pixels included in the peak
#define CAT 11      // peak category [0..9]

#define FITI 12     // fit iterations
#define FRAME 13    // frame the peak first occured in
#define TLEN 14     // track length
#define LINK 15     // id of the next molecule in the track

#define ZO 16       // original z coordinate
#define Z 17        // drift corrected z coordinate

