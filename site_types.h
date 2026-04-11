enum SiteType {
    TOP,
    HOLLOW,
    BRIDGE,
    // Other existing site types
};

class SiteTypeManager {
public:
    static SiteType parseSiteType(const std::string& input) {
        // Implement parsing logic for input reactions like "A + t -> A(t)"
        if (input == "TOP") return SiteType::TOP;
        else if (input == "HOLLOW") return SiteType::HOLLOW;
        else if (input == "BRIDGE") return SiteType::BRIDGE;
        // Additional parsing logic...
        throw std::invalid_argument("Invalid site type");
    }

    // Add properties management methods as needed
};
