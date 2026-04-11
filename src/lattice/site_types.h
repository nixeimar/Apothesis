# Site Types

This header file defines the enumeration and structures for different site types in a lattice structure.

## Definitions

### Site Types Enumeration
```c
typedef enum {
    TOP,
    HOLLOW,
    BRIDGE
} SiteType;
```

### Site Structure
```c
typedef struct {
    SiteType type;
    double properties[3]; // to hold site-specific properties
    int coordination_number;
} Site;
```

// Function prototypes and other declarations can be added here as needed.
