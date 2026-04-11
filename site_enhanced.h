#ifndef SITE_ENHANCED_H
#define SITE_ENHANCED_H

enum SiteType { TOP, HOLLOW, BRIDGE };

class Site {
public:
    Site(SiteType type);
    void displayType();
private:
    SiteType type_;
};

#endif // SITE_ENHANCED_H
