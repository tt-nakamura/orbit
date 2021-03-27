#ifndef __Observatory_h__
#define __Observatory_h__

struct Observatory {
    const char *code;
    double longitude;
    double latitude;
    static Observatory *lookup(const char*);
};

#endif // __Observatory_h__