// rx_phys.cpp - MIT License
// Roman Kobylinsky 2021
//
// This is a simple 2d physics engine written in a single weekend.
// Hope that someone will find this project useful.
//
// Project depends on imgui (docking branch) and sokol libraries, and requires C++20 to compile.
//
// Main recipe of this physics engine is:
//      - Advance bodies along their trajectories
//      - Find all intersecting pairs
//      - For all intersecting pairs
//          - Push them out of each other
//          - Update their velocity based on mass and restitution
//      - ???
//      - Repeat
//
// For finding a intersection between bodies, I chose general GJK algorithm,
// which is a bit complicated but allows for easy addition of new types of bodies.
// But result of GJK algorithm does not yeild enough information, that's why
// there's a second part: An EPA algorithm. Which helps us find penetration depth and normal.
// As for updaing velocities, it is basically seventh grade physics class,
// though you need to be careful with numerical errors.
//
// Further work:
//      - Cache everything that can be cached
//      - Add special case for spheres in GJK
//          - Sphere collision can go without EPA
//      - Rotation (this one actually seems very simple and cool to do)
//      - Friction (other than rotation seems pretty important to have)
//      - Kinematic bodies (currently objects with infinite mass count as that)
//      - Use BVH acceleration scructure
//      - SOA all the things
//      - Spring joints
//      - GJK and EPA introspection and visualization
//      - 3D (3 main annoyances)
//          - GJK simplex case
//          - EPA for polyhedra
//          - Angular momentum, seems scary (probably not)
//          - https://www.youtube.com/watch?v=1n-HMSCDYtM :)
//
// Index of this file:
// [SECTION] Header mess
// [SECTION] Common utilities
// [SECTION] Bodies declaration
// [SECTION] Body Implementation
// [SECTION] Circle Implementation
// [SECTION] Box Implementation
// [SECTION] Polygon Implementation
// [SECTION] Gilbert Johnson Keerthi Stelly Muratori algorithm :)
// [SECTION] Expanding polytope algorithm.
// [SECTION] Collision resolution
// [SECTION] World declaration
// [SECTION] World implementation
// [SECTION] Entry point, more or less...

//==-----------------------------------------------------------------------==//
// [SECTION] Header mess
//==-----------------------------------------------------------------------==//

#include <imgui.h>
#include <imgui_internal.h>
#include <memory>
#include <vector>

#include <sokol_app.h>
#include <sokol_gfx.h>
#include <sokol_glue.h>
#include <sokol_time.h>
#include <sokol_log.h>

#include <sokol.h>
#define SOKOL_IMGUI_IMPL
#include <util/sokol_imgui.h>

//==-----------------------------------------------------------------------==//
// [SECTION] Common utilities
//==-----------------------------------------------------------------------==//

#define MIN(a, b) ((a) < (b)) ? (a) : (b)
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define CLAMP(v, min, max) MAX(min, MIN(max, v))

inline ImVec2 Normalize(ImVec2 vector) {
    ImVec2 result = vector / ImSqrt(ImLengthSqr(vector));
    return result;
}

bool SameDirection(ImVec2 a, ImVec2 b) {
    bool result = ImDot(a, b) > 0.0f;
    return result;
}

float Cross(ImVec2 a, ImVec2 b) {
    float result = a.x * b.y - a.y * b.x;
    return result;
}

inline bool PointOnTheRight(ImVec2 a, ImVec2 b, ImVec2 point) {
    float dotProduct = ((b.y - a.y) * (point.x - a.x)) + ((a.x - b.x) * (point.y - a.y));
    return (dotProduct > 0.0f);
}

inline ImVec2 Reflect(ImVec2 a, ImVec2 normal) {
    float dn = 2 * ImDot(a, normal);
    return a - normal * dn;
}

//==-----------------------------------------------------------------------==//
// [SECTION] Bodies declaration
//==-----------------------------------------------------------------------==//

/**
 * @brief Abstract representation of different bodies and their properties.
 */
class Body {
  public:
    /**
     * @brief Computes axis align boudning box using support function.
     */
    ImRect ComputeAABB();

    /**
     * @brief Returns farthest point on the rim of the body in given direction.
     * @param normal In which direction to search for, don't have to be normalized.
     */
    virtual ImVec2 Support(ImVec2 normal) = 0;
    /**
     * @brief Render the body within provided imgui context.
     */
    virtual void Render(ImDrawList *drawList, ImVec2 origin, ImU32 color,
                        float pixelsPerMeter) = 0;
    virtual const char *Name() = 0; ///< Returns name of the object.

    /**
     * @brief Advance this particular object in time, applying gravity and moving along it's velocity.
     * @param deltaTime Amount of time in seconds.
     */
    void Advance(float deltaTime, ImVec2 gravity);

    void SetPosition(ImVec2 position); ///< Set Position in world space.
    ImVec2 GetPosition() const; ///< Returns current position in world space.

    void SetVelocity(ImVec2 velocity); ///< Set velocity in m/s.
    ImVec2 GetVelocity() const;        ///< Returns current velocity in m/s.

    /**
     * @brief Set mass of the object in kg.
     * Zero is interpreted as infinite mass.
     */
    void SetMass(float mass);
    /**
     * @brief Returns mass of the object in kg.
     * Zero is interpreted as infinite mass.
     */
    float GetMass() const;

    void SetRestitution(float restitution); ///< Set restitution of the object.
    float GetRestitution() const;           ///< Returns current restitution.

  protected:
    /**
     * @brief Construct a new Body object.
     *
     * @param mass Mass of the object in kilograms.
     * @param restitution Amount of bounciness of the object in range from 0 to 1.
     * @param position Position in meters in world space.
     */
    Body(float mass, float restitution, ImVec2 position);

    float m_mass = 1.0f;        ///< Mass of the object in killograms.
    float m_restitution = 0.5f; ///< Bounciness of the object in range [0, 1].
    ImVec2 m_position = {};     ///< Position in world space in meters.
    ImVec2 m_velocity = {};     ///< Current velocity in meters per second.
};

class Circle : public Body {
  public:
    /**
     * @brief Construct a new Circle object.
     *
     * @param mass Mass of the object in kilograms.
     * @param restitution Amount of bounciness of the object in range from 0 to 1.
     * @param position Position in meters in world space.
     * @param radius Radius of the circle in meters.
     */
    Circle(float mass, float restitution, ImVec2 position, float radius);

    /**
     * @brief Returns farthest point on the rim of the body in given direction.
     * @param normal In which direction to search for, don't have to be normalized.
     */
    ImVec2 Support(ImVec2 normal) override;
    /**
     * @brief Render the body within provided imgui context.
     */
    void Render(ImDrawList *drawList, ImVec2 origin, ImU32 color, float pixelsPerMeter) override;
    const char *Name() override; ///< Returns name of the object.

  private:
    float m_radius; ///< Radius of the circle in meters.
};

class Box : public Body {
  public:
    /**
     * @brief Construct a new Box object.
     *
     * @param mass Mass of the object in kilograms.
     * @param restitution Amount of bounciness of the object in range from 0 to 1.
     * @param position Position in meters in world space.
     * @param radius Radius of the box in meters. Actual dimensions of the box would be twice as big!
     */
    Box(float mass, float restitution, ImVec2 position, ImVec2 radius);

    /**
     * @brief Returns farthest point on the rim of the body in given direction.
     * @param normal In which direction to search for, don't have to be normalized.
     */
    ImVec2 Support(ImVec2 normal) override;
    /**
     * @brief Render the body within provided imgui context.
     */
    void Render(ImDrawList *drawList, ImVec2 origin, ImU32 color, float pixelsPerMeter) override;
    const char *Name() override; ///< Returns name of the object.

  private:
    ImVec2 m_radius; ///< Radius of the box in meters.
};

class Poly : public Body {
  public:
    /**
     * @brief Construct a new Polygon object.
     *
     * @param mass Mass of the object in kilograms.
     * @param restitution Amount of bounciness of the object in range from 0 to 1.
     * @param position Position in meters in world space.
     * @param vertices Vertices of the polygons in local space and in meters.
     */
    Poly(float mass, float restitution, ImVec2 position, std::vector<ImVec2> vertices);

    /**
     * @brief Returns farthest point on the rim of the body in given direction.
     * @param normal In which direction to search for, don't have to be normalized.
     */
    ImVec2 Support(ImVec2 normal) override;
    /**
     * @brief Render the body within provided imgui context.
     */
    void Render(ImDrawList *drawList, ImVec2 origin, ImU32 color, float pixelsPerMeter) override;
    const char *Name() override; ///< Returns name of the object.

  private:
    std::vector<ImVec2> m_vertices; ///< Vertices in local space.
};

//==-----------------------------------------------------------------------==//
// [SECTION] Body Implementation
//==-----------------------------------------------------------------------==//

inline Body::Body(float mass, float restitution, ImVec2 position)
    : m_mass(mass), m_restitution(restitution), m_position(position) {}

inline ImRect Body::ComputeAABB() {
    ImRect result;
    result.Min.x = Support(ImVec2(-1.0f, 0.0f)).x;
    result.Min.y = Support(ImVec2(0.0f, -1.0f)).y;
    result.Max.x = Support(ImVec2(1.0f, 0.0f)).x;
    result.Max.y = Support(ImVec2(0.0f, 1.0f)).y;
    return result;
};

inline void Body::Advance(float deltaTime, ImVec2 gravity) {
    // TODO(kobylinsky): Split 'infinite mass' and 'stationary' into different properties.
    if (m_mass != 0.0f) {
        m_position += m_velocity * deltaTime + gravity * (deltaTime * deltaTime) * 0.5f;
        m_velocity += gravity * deltaTime;
    }
}

inline void Body::SetPosition(ImVec2 position) { m_position = position; }
inline ImVec2 Body::GetPosition() const { return m_position; }

inline void Body::SetVelocity(ImVec2 velocity) { m_velocity = velocity; }
inline ImVec2 Body::GetVelocity() const { return m_velocity; }

inline void Body::SetMass(float mass) { m_mass = mass; }
inline float Body::GetMass() const { return m_mass; }

inline void Body::SetRestitution(float restitution) {
    m_restitution = restitution;
}
inline float Body::GetRestitution() const { return m_restitution; }

//==-----------------------------------------------------------------------==//
// [SECTION] Circle Implementation
//==-----------------------------------------------------------------------==//

inline Circle::Circle(float mass, float restitution, ImVec2 position, float radius)
    : Body(mass, restitution, position), m_radius(radius){};

inline ImVec2 Circle::Support(ImVec2 normal) {
    return m_position + Normalize(normal) * m_radius;
}

inline void Circle::Render(ImDrawList *drawList, ImVec2 origin, ImU32 color, float pixelsPerMeter) {
    drawList->AddCircle(origin + m_position * pixelsPerMeter, m_radius * pixelsPerMeter, color);
}

inline const char *Circle::Name() { return "Circle"; }

//==-----------------------------------------------------------------------==//
// [SECTION] Box Implementation
//==-----------------------------------------------------------------------==//

inline Box::Box(float mass, float restitution, ImVec2 position, ImVec2 radius)
    : Body(mass, restitution, position), m_radius(radius) {};

inline ImVec2 Box::Support(ImVec2 normal) {
    ImVec2 result;
    result.x = (normal.x > 0.0f) ? m_radius.x : -m_radius.x;
    result.y = (normal.y > 0.0f) ? m_radius.y : -m_radius.y;
    return m_position + result;
}

inline void Box::Render(ImDrawList *drawList, ImVec2 origin, ImU32 color, float pixelsPerMeter) {
    drawList->AddRect(origin + (m_position - m_radius) * pixelsPerMeter,
                      origin + (m_position + m_radius) * pixelsPerMeter, color);
}

inline const char *Box::Name() { return "Box"; }

//==-----------------------------------------------------------------------==//
// [SECTION] Polygon Implementation
//==-----------------------------------------------------------------------==//

inline Poly::Poly(float mass, float restitution, ImVec2 position,
                  std::vector<ImVec2> vertices)
    : Body(mass, restitution, position), m_vertices(vertices){};

inline ImVec2 Poly::Support(ImVec2 normal) {
    ImVec2 result = {};
    float maxDistance = -FLT_MAX;
    for (auto &&vertex : m_vertices) {
        auto candidate = m_position + vertex;
        float distance = ImDot(candidate, normal);
        if (distance > maxDistance) {
            maxDistance = distance;
            result = candidate;
        }
    }
    return result;
}

inline void Poly::Render(ImDrawList *drawList, ImVec2 origin, ImU32 color, float pixelsPerMeter) {
    std::vector<ImVec2> points;
    for (auto &&vertex : m_vertices) {
        points.emplace_back(origin + (m_position + vertex) * pixelsPerMeter);
    }
    if (!points.empty()) {
        points.emplace_back(points[0]);
    }
    drawList->AddPolyline(points.data(), (int)points.size(), color, 0, 1.0f);
}
inline const char *Poly::Name() { return "Poly"; }

//==-----------------------------------------------------------------------==//
// [SECTION] Gilbert Johnson Keerthi Stelly Muratori algorithm :)
//==-----------------------------------------------------------------------==//

struct Simplex {
    size_t maxIterations = 20;   ///< How many iteration to try before bailing.
    size_t currentIteration = 0; ///< Current iteration, updated by @ref GJK

    bool hit = false;       ///< Whether hit were dectected.
    size_t vertexCount = 0; ///< Number of vertices in the simplex.

    /**
     * @brief Vertices of the simplex, actual amount is stored in @ref
     * Simplex::vertexCount
     */
    struct Vertex {
        // TODO(kobylinsky): Move this structure out of @ref Simplex.
        ImVec2 a = {}; ///< Vertex on the rim of the first body.
        ImVec2 b = {}; ///< Vertex on the rim of the second body.
        ImVec2 p = {}; ///< Vertex on the minkowski difference of two bodies.
    } vertices[3] = {};

    /**
     * @brief Intermediate values used to compute next query vector.
     */
    float barycenter[3] = {};
    float distanceSquared = FLT_MAX; ///< Squared distance between two bodies.
};

/**
 * @brief Structure that is used to communicate with @ref GJK and @ref EPA
 * algorithms. @ref GJK and @ref EPA after each iteration will fill out which
 * direction they want to sample, and user has to provide according vertices
 * before the next iteration.
 */
struct Support {
    // Input
    ImVec2 aVertex = {}; ///< Farthest vertex along @ref aDirection.
    ImVec2 bVertex = {}; ///< Farthest vertex along @ref bDirection.

    // Output
    ImVec2 aDirection = {}; ///< Direction in which to sample A body.
    ImVec2 bDirection = {}; ///< Direction in which to sample B body.
};

/**
 * @brief GJK simplex algorithm. Look below on how to use it.
 * Reference https://caseymuratori.com/blog_0003
 *
 * @param simplex Intermediate state of the algorithm.
 * @param support Stucture used to communicate with the user.
 * @return true Algorithm is still in progress, update @ref Support sturcture and call it again.
 * @return false Algorithm has done it's work, now by examining @ref Simplex you could find whether hit occured, or get distance between objects otherwise.
 *
 * Here's how you use this function:
 * @code
 * // Start by providint two initial points.
 * Support support = {.aVertex = a.Support(ImVec2(-1, 0)),
 *                    .bVertex = b.Support(ImVec2( 1, 1))};
 * // Initialize the simplex structure
 * Simplex simplex = {
 *     .maxIterations = 20, // Here's how you can limit number of iterations.
 * };
 * // Next iterate unitl it returns true.
 * while (GJK(&simplex, &support)) {
 *     // Provide requested verticies.
 *     support.aVertex = a.Support(support.aDirection);
 *     support.bVertex = b.Support(support.bDirection);
 * }
 * // Now we can examine the simplex.
 * if (simplex.hit) {
 *     // Yay, we got collision.
 * } else {
 *     // Objects did not collide.
 *     // But here's the distance between them:
 *     float distance = sqrtf(simplex.distanceSquared);
 * }
 * @endcode
 */
bool GJK(Simplex *simplex, Support *support) {
    assert(simplex && "Simplex is not supplied.");
    assert(support && "Support is not supplied.");
    assert((simplex->vertexCount <= 3) && "Too many vertices!");

    if (!simplex || !support) {
        return false;
    }

    if ((simplex->maxIterations > 0) &&
        (simplex->currentIteration++ >= simplex->maxIterations)) {
        return false;
    }

    // Copy vertex from 'Support' into simplex.
    simplex->vertices[simplex->vertexCount] = {
        .a = support->aVertex,
        .b = support->bVertex,
        .p = support->bVertex - support->aVertex,
    };
    simplex->barycenter[simplex->vertexCount++] = 1.0f;

    ImVec2 searchDirection = -simplex->vertices[0].p;

    switch (simplex->vertexCount) {
        case 1:
            break;
        case 2: {
            auto a = simplex->vertices[0].p;
            auto b = simplex->vertices[1].p;

            auto ab = a - b;
            auto ba = b - a;

            float u = ImDot(b, ba);
            float v = ImDot(a, ab);

            if (v <= 0.0f) {
                simplex->barycenter[0] = 1.0f;
                simplex->vertexCount = 1;
                break;
            }
            if (u <= 0.0f) {
                simplex->vertices[0] = simplex->vertices[1];
                simplex->barycenter[0] = 1.0f;
                simplex->vertexCount = 1;
                break;
            }
            simplex->barycenter[0] = u;
            simplex->barycenter[1] = v;
            simplex->vertexCount = 2;
        } break;
        case 3: {
            auto a = simplex->vertices[0].p;
            auto b = simplex->vertices[1].p;
            auto c = simplex->vertices[2].p;

            auto ab = b - a;
            auto ca = c - a;
            auto bc = c - b;

            float v_ab = -ImDot(a, ab);
            float u_ab =  ImDot(b, ab);

            float u_ca = -ImDot(a, ca);
            float v_ca =  ImDot(c, ca);

            float v_bc = -ImDot(b, bc);
            float u_bc =  ImDot(c, bc);

            if (v_ab <= 0.0f && u_ca <= 0.0f) {
                simplex->barycenter[0] = 1.0f;
                simplex->vertexCount = 1;
                break;
            }

            if (u_ab <= 0.0f && v_bc <= 0.0f) {
                simplex->barycenter[1] = 1.0f;
                simplex->vertexCount = 1;

                simplex->vertices[0] = simplex->vertices[1];
                simplex->barycenter[0] = simplex->barycenter[1];
                break;
            }

            if (v_ca <= 0.0f && u_bc <= 0.0f) {
                simplex->barycenter[2] = 1.0f;
                simplex->vertexCount = 1;

                simplex->vertices[0] = simplex->vertices[2];
                simplex->barycenter[0] = simplex->barycenter[2];
                break;
            }

            float n = Cross(ab, ca);

            float u_abc = n * Cross(b, c);
            float v_abc = n * Cross(c, a);
            float w_abc = n * Cross(a, b);

            if (u_ab > 0.0f && v_ab > 0.0f && w_abc <= 0.0f) {
                float inv_ab = 1.0f / (u_ab + v_ab);
                simplex->barycenter[0] = u_ab * inv_ab;
                simplex->barycenter[1] = v_ab * inv_ab;
                simplex->vertexCount = 2;
                break;
            }

            if (v_ca > 0.0f && u_ca > 0.0f && v_abc <= 0.0f) {
                float inv_ac = 1.0f / (v_ca + u_ca);
                simplex->barycenter[0] = v_ca * inv_ac;
                simplex->barycenter[2] = u_ca * inv_ac;
                simplex->vertexCount = 2;

                simplex->vertices[1] = simplex->vertices[2];
                simplex->barycenter[1] = simplex->barycenter[2];

                break;
            }

            if (u_bc > 0.0f && v_bc > 0.0f && u_abc <= 0.0f) {
                float inv_bc = 1.0f / (u_bc + v_bc);
                simplex->barycenter[1] = u_bc * inv_bc;
                simplex->barycenter[2] = v_bc * inv_bc;
                simplex->vertexCount = 2;

                simplex->vertices[0] = simplex->vertices[2];
                simplex->barycenter[0] = simplex->barycenter[2];
                break;
            }

            float inv_abc = 1.0f / (u_abc + v_abc + w_abc);
            simplex->barycenter[0] = u_abc * inv_abc;
            simplex->barycenter[1] = v_abc * inv_abc;
            simplex->barycenter[2] = w_abc * inv_abc;
            simplex->vertexCount = 3;

            simplex->hit = 1;
            return false;
        } break;
    }

    ImVec2 point;
    float denom = 0.0f;

    for (size_t i = 0; i < simplex->vertexCount; ++i) {
        denom += simplex->barycenter[i];
    }
    denom = 1.0f / denom;

    switch (simplex->vertexCount) {
    case 1:
        point = simplex->vertices[0].p;
        break;
    case 2: {
        ImVec2 a = simplex->vertices[0].p * (denom * simplex->barycenter[0]);
        ImVec2 b = simplex->vertices[1].p * (denom * simplex->barycenter[1]);
        point = a + b;
    } break;
    }

    float d2 = ImLengthSqr(point);
    if (d2 >= simplex->distanceSquared) {
        return false;
    }
    simplex->distanceSquared = d2;

    ImVec2 d = {};
    switch (simplex->vertexCount) {
    default:
        assert(false);
        break;
    case 1: {
        d = -simplex->vertices[0].p;
    } break;
    case 2: {
        auto a = simplex->vertices[0].p;
        auto b = simplex->vertices[1].p;
        auto ab = b - a;

        if (PointOnTheRight(a, b, ImVec2())) {
            d = ImVec2(ab.y, -ab.x);
        } else {
            d = ImVec2(-ab.y, ab.x);
        }
    } break;
    }

    if (ImLengthSqr(d) < FLT_EPSILON) {
        return false;
    }

    support->aDirection = -d;
    support->bDirection =  d;

    return true;
}

//==-----------------------------------------------------------------------==//
// [SECTION] Expanding polytope algorithm.
//==-----------------------------------------------------------------------==//

struct Polytope {
    /**
     * @brief Constructs initial polytope from simplex resulting from GJK
     * algorithm.
     */
    Polytope(Simplex *simplex, size_t maxIterations = 20)
        : maxIterations(maxIterations) {
        assert(simplex->hit && "Simplex should include the origin!");
        for (size_t vertexIndex = 0; vertexIndex < simplex->vertexCount;
             ++vertexIndex) {
            vertices.emplace_back(simplex->vertices[vertexIndex]);
        }
    }

    size_t maxIterations = 20;   ///< How many iteration to try before bailing.
    size_t currentIteration = 0; ///< Current iteration, updated by @ref GJK

    std::vector<Simplex::Vertex> vertices; ///< Vertices of the polytope.

    float distance = INFINITY; ///< Distance from polytope to the origin.
    size_t index = 0;          ///< Index of the closest vertex to the origin.
    ImVec2 normal = {};        ///< Normal of the closest edge to the origin.

    bool hit = false;        ///< Wheter result was found.
    ImVec2 penetration = {}; /// < Penetration vector.

    /**
     * @brief First iteration does not know which direction it wants to sample.
     * So here's a little latch to ignore first iteration when adding vertex to
     * the polytope.
     */
    bool firstIteration = true;
};

/**
 * @brief Expanding polytope algorithm - Grows initial polytope in the direction
 * of the origin. Reference https://dyn4j.org/2010/05/epa-expanding-polytope-algorithm/
 *
 * @param polytope Intermediate state of the algorithm.
 * @param support Stucture used to communicate with the user.
 * @return true Algorithm is still in progress, update @ref Support sturcture and call it again.
 * @return false Algorithm has done it's work, now by examining @ref Polytope you could find collision normal and penetration depth.
 *
 * Here's how you use this function:
 * @code
 * // Construct initial polytope from the simplex resulting from GJK.
 * // Note that simplex should have 'hit' set to true.
 * Polytope polytope(&simplex, 20);
 * // Initialize the support structure
 * Support support = {};
 * while (EPA(&polytope, &support)) {
 *     // Provide requested support vertices.
 *     support.aVertex = a.Support(support.aDirection);
 *     support.bVertex = b.Support(support.bDirection);
 * }
 * if (polytope.hit) {
 *     // Now we have enough information about the collision.
 *     ResolveNarrowCollision(a, b, polytope.normal);
 * }
 * @endcode
 */
inline bool EPA(Polytope *polytope, Support *support) {
    assert(polytope && "Simplex is not supplied.");
    assert(support && "Support is not supplied.");
    assert((polytope->vertices.size() >= 3) && "Too few vertices!");

    auto &p = *polytope;

    if ((polytope->maxIterations > 0) &&
        (polytope->currentIteration++ >= polytope->maxIterations)) {
        return false;
    }

    if (p.firstIteration) {
        p.firstIteration = false;
    } else {
        auto supportPoint = support->bVertex - support->aVertex;
        auto distance = ImDot(p.normal, supportPoint);
        if (abs(distance - p.distance) > 0.001) {
            p.distance = INFINITY;
            p.vertices.insert(p.vertices.begin() + p.index,
                              Simplex::Vertex{
                                  .a = support->aVertex,
                                  .b = support->bVertex,
                                  .p = supportPoint,
                              });
        }
        if (p.distance != INFINITY) {
            p.hit = true;
            p.penetration = p.normal * (p.distance + FLT_EPSILON);
            return false;
        }
    }

    for (size_t vertexIndex = 0; vertexIndex < p.vertices.size();
         vertexIndex++) {
        size_t nextVertexIndex = (vertexIndex + 1) % p.vertices.size();

        auto a = p.vertices[vertexIndex].p;
        auto b = p.vertices[nextVertexIndex].p;

        auto ab = b - a;

        auto normal = Normalize(ImVec2(ab.y, -ab.x));
        auto distance = ImDot(normal, a);

        if (distance < 0) {
            distance *= -1;
            normal *= -1;
        }

        if (distance < p.distance) {
            p.distance = distance;
            p.normal = normal;
            p.index = nextVertexIndex;
        }
    }

    support->aDirection = -p.normal;
    support->bDirection =  p.normal;

    return true;
}

//==-----------------------------------------------------------------------==//
// [SECTION] Collision resolution
//==-----------------------------------------------------------------------==//

/**
 * @brief Update speeds of the objects, assuming they are collided.
 * @param normal Normal of the contact.
 */
inline void ResolveNarrowCollision(Body &a, Body &b, ImVec2 normal) {
    ImVec2 relativeVelocity = b.GetVelocity() - a.GetVelocity();
    float velocityAlongNormal = ImDot(relativeVelocity, normal);

    if (velocityAlongNormal > 0)
        return; // Bail if objects moving away from each other.

    float restitution = MIN(a.GetRestitution(), b.GetRestitution());

    // TODO(kobylinsky): Cache inverse mass in @ref Body object.
    float totalMass = a.GetMass() + b.GetMass();
    float aRatio = a.GetMass() / totalMass;
    float bRatio = b.GetMass() / totalMass;

    float aInvMass = 0.0f;
    float bInvMass = 0.0f;
    if (a.GetMass() != 0.0f) {
        aInvMass = 1.0f / a.GetMass();
    }
    if (b.GetMass() != 0.0f) {
        bInvMass = 1.0f / b.GetMass();
    }

    float power =
        ((1 + restitution) * -velocityAlongNormal) / (aInvMass + bInvMass);
    ImVec2 impulse = normal * power;
    a.SetVelocity(a.GetVelocity() - impulse * aRatio);
    b.SetVelocity(b.GetVelocity() + impulse * bRatio);
}

/**
 * @brief Figure out whether two bodies are colliding and update their velocity
 * if they are.
 */
inline void ResolveCollision(Body &a, Body &b) {
    // If two bodies have an infinite mass it is pointless to resolve it.
    // They should pass right through each other.
    if (a.GetMass() == 0 && b.GetMass() == 0) {
        return;
    }

    // If AABBs of the objects do not overlap, there's no collision.
    // TODO(kobylinsky): Move AABB check somewhere, where colliding pairs will
    // be searched for.
    if (!a.ComputeAABB().Overlaps(b.ComputeAABB())) {
        return;
    }

    // @ref Support structure is used to communicate with the algorithm
    // while leaving control flow to the user.
    Support support = {.aVertex = a.Support(ImVec2(-1, 0)),
                       .bVertex = b.Support(ImVec2(1, 1))};

    Simplex simplex = {};
    // First we use GJK to figure out if there's an intersection.
    while (GJK(&simplex, &support)) {
        // Provide requested support vertices.
        support.aVertex = a.Support(support.aDirection);
        support.bVertex = b.Support(support.bDirection);
    }
    if (simplex.hit) {
        // As a result we get a simplex enclosing the origin
        // of the minkowski difference. In order to figure out actual
        // penetration vector we need to find an edge of the minkowski
        // difference that is closest to the origin. We do that by leveraging
        // Expanding Polytope algorithm, that will grow our find simplex in the
        // direction of the origin.

        // First convert newfound simplex into a polytope.
        Polytope polytope(&simplex);
        Support support = {};
        while (EPA(&polytope, &support)) {
            // Provide requested support vertices.
            support.aVertex = a.Support(support.aDirection);
            support.bVertex = b.Support(support.bDirection);
        }
        if (polytope.hit) {
            // Now we have enough information about the collision to
            // update velocities of our bodies...
            ResolveNarrowCollision(a, b, -polytope.normal);
            // ...and push them out of each other.
            // NOTE(kobylinsky): We don't push out static objects.
            // TODO(kobylinsky): Split 'infinite mass' and 'stationary' into
            // different properties.
            if (a.GetMass() != 0.0f)
                a.SetPosition(a.GetPosition() + polytope.penetration * 0.5f);
            if (b.GetMass() != 0.0f)
                b.SetPosition(b.GetPosition() - polytope.penetration * 0.5f);
        }
    } else {
        // Even if objects do not collide, we can still extract some useful
        // information out of the simplex. Like the distance.
        // TODO(kobylinsky): Do the distance function between bodies.
    }
}

//==-----------------------------------------------------------------------==//
// [SECTION] World declaration
//==-----------------------------------------------------------------------==//

class World {
  public:
    World() = default;

    template <typename BodyType, typename... Args>
    std::shared_ptr<BodyType> AddBody(Args &&...args) {
        static_assert(std::is_base_of<Body, BodyType>::value, "All physics bodies must be inherited from Body!");
        auto body = std::make_shared<BodyType>(std::forward<Args>(args)...);
        m_bodies.emplace_back(body);
        return body;
    }

    void Update(float deltaTime);
    void WindowConfig(const char *label);
    void WindowViewport(const char *label, ImVec2 &scrolling);

  private:
    void ContextMenu(ImVec2 mouseInCanvas);
    void DrawGrid(ImDrawList *drawList, ImRect canvasBounds, ImVec2 scrolling);

    std::vector<std::shared_ptr<Body>> m_bodies;
    float m_timeScale = 1.0f;
    ImVec2 m_gravity = {0.0f, 9.81f};
    float m_pixelsPerMeter = 64.0f;
    bool m_drawBoundingBoxes = false;
    bool m_drawVelocities = false;
    bool m_update = true;
    float m_largestTimeStep = 0.01f;
};

//==-----------------------------------------------------------------------==//
// [SECTION] World implementation
//==-----------------------------------------------------------------------==//

inline void World::Update(float deltaTime) {
    if (!m_update) {
        return;
    }

    // Find out how much time we are actually simulating.
    float time = deltaTime * m_timeScale;

    while (time > 0) {
        // Step should be no larger than provided, so we split
        // our time into chunks and simulate each chunk separately.
        float step = MIN(time, m_largestTimeStep);

        // Advance gravity, velocity and translation of all bodies.
        for (auto &&body : m_bodies) {
            body->Advance(step, m_gravity);
        }

        // For each pair of objects try to resolve collision.
        for (size_t i = 0; i < m_bodies.size(); ++i) {
            for (size_t j = i + 1; j < m_bodies.size(); ++j) {
                ResolveCollision(*m_bodies[i], *m_bodies[j]);
            }
        }
        time -= step;
    }
}

inline void World::WindowConfig(const char *label) {
    using namespace ImGui;

    SetNextWindowSize(ImVec2(660, 370), ImGuiCond_Appearing);
    if (Begin(label)) {
        SliderFloat("Time Scale", &m_timeScale, 0.0f, 5.0f);
        SliderFloat("Largest Time Step", &m_largestTimeStep, 0.001f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
        DragFloat2("Gravity", &m_gravity.x, -20.0f, 20.0f);
        SliderFloat("Pixels per meter", &m_pixelsPerMeter, 10.0f, 256.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
        Checkbox("Bounding Boxes", &m_drawBoundingBoxes);
        Checkbox("Velocity Vectors", &m_drawVelocities);
        Checkbox("Update", &m_update);

        if (Button("Reset")) {
            m_bodies.clear();

            AddBody<Box>(0.0f, 0.5f, ImVec2(0, -10.0f), ImVec2(10.0f, 0.25f));
            AddBody<Box>(0.0f, 0.5f, ImVec2(-10.0f, 0), ImVec2(0.25f, 10.0f));
            AddBody<Box>(0.0f, 0.5f, ImVec2(10.0f, 0), ImVec2(0.25f, 10.0f));
            AddBody<Box>(0.0f, 0.5f, ImVec2(0, 10.0f), ImVec2(10.0f, 0.25f));
        }
        SameLine();
        Text("Simulating %d bodies.", m_bodies.size());

        if (Button("Earth")) {
            m_gravity = ImVec2(0.0f, 9.807f);
        }
        SameLine();
        if (Button("Mars")) {
            m_gravity = ImVec2(0.0f, 3.721f);
        }
        SameLine();
        if (Button("Moon")) {
            m_gravity = ImVec2(0.0f, 1.62f);
        }
        SameLine();
        if (Button("Zero")) {
            m_gravity = ImVec2(0.0f, 0.0f);
        }

        if (Button("Real-time")) {
            m_timeScale = 1.0f;
        }
        SameLine();
        if (Button("Slow")) {
            m_timeScale = 0.001f;
        }
        SameLine();
        if (Button("Stop")) {
            m_timeScale = 0.0f;
        }
        SameLine();
        if (Button("Step")) {
            m_timeScale = 1.0f;
            Update(0.016f);
            m_timeScale = 0.0f;
        }

        Separator();

        TextWrapped(
            "This is a simple physics engine, which was done within a single weekend, which I enjoyed a whole lot!\n"
            "This is an inspector window which exposes various parameters of the system.\n"
            "There's also a viewport window:"
        );
        BulletText("Use Left mouse to pick and throw physics bodies");
        BulletText("Click Right mouse to open object spawn menu");
        BulletText("By holding right mouse button you can move the camera");
        TextWrapped("Check the source code too! It is even somewhat documented!");
    }
    End();
}

inline void World::WindowViewport(const char *label, ImVec2 &scrolling) {
    using namespace ImGui;

    if (Begin(label)) {
        ImRect canvasBounds;
        canvasBounds.Min = GetCursorScreenPos();
        canvasBounds.Max = canvasBounds.Min + ImMax(GetContentRegionAvail(), ImVec2(50, 50));

        ImGuiIO &io = GetIO();
        ImDrawList *drawList = GetWindowDrawList();
        drawList->AddRectFilled(canvasBounds.Min, canvasBounds.Max, IM_COL32(50, 50, 50, 255));
        drawList->AddRect(canvasBounds.Min, canvasBounds.Max, IM_COL32(255, 255, 255, 255));

        InvisibleButton("canvas", canvasBounds.GetSize(), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
        const bool isHovered = IsItemHovered();
        const bool isActive = IsItemActive();
        const ImVec2 origin(canvasBounds.Min.x + scrolling.x, canvasBounds.Min.y + scrolling.y);
        const ImVec2 mouseInCanvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

        if (isActive && IsMouseDragging(ImGuiMouseButton_Right, -1.0f)) {
            scrolling.x += io.MouseDelta.x;
            scrolling.y += io.MouseDelta.y;
        }
        ImVec2 dragDelta = GetMouseDragDelta(ImGuiMouseButton_Right);
        m_pixelsPerMeter =
            CLAMP(m_pixelsPerMeter + GetIO().MouseWheel, 10.0f, 256.0f);

        if (IsMouseReleased(ImGuiMouseButton_Right) && dragDelta.x == 0.0f && dragDelta.y == 0.0f) {
            OpenPopupOnItemClick("context");
        }
        if (BeginPopup("context")) {
            ContextMenu(mouseInCanvas);
            EndPopup();
        }

        static std::weak_ptr<Body> selectedWeak;
        std::shared_ptr<Body> selected = selectedWeak.lock();

        if (IsMouseClicked(ImGuiMouseButton_Left)) {
            for (auto body : m_bodies) {
                if (body->ComputeAABB().Contains(mouseInCanvas / m_pixelsPerMeter)) {
                    selected = body;
                    selectedWeak = selected;
                    break;
                }
            }
        }

        if (isActive && IsMouseDown(ImGuiMouseButton_Left) && selected) {
            selected->SetPosition(mouseInCanvas / m_pixelsPerMeter);
            if (m_timeScale != 0.0f) {
                selected->SetVelocity((GetIO().MouseDelta / (GetIO().DeltaTime * m_timeScale)) / m_pixelsPerMeter);
            } else {
                selected->SetVelocity(ImVec2());
            }
        }

        if (!isActive || IsMouseReleased(ImGuiMouseButton_Left)) {
            selected.reset();
            selectedWeak.reset();
        }

        drawList->PushClipRect(canvasBounds.Min, canvasBounds.Max, true);
        DrawGrid(drawList, canvasBounds, scrolling);
        { // Draw bodies to the canvas.
            for (auto body : m_bodies) {
                body->Render(drawList, origin, (body == selected) ? IM_COL32(255, 0, 0, 255) : IM_COL32(0, 255, 255, 255), m_pixelsPerMeter);
                if (m_drawBoundingBoxes) {
                    auto aabb = body->ComputeAABB();
                    drawList->AddRect(origin + aabb.Min * m_pixelsPerMeter, origin + aabb.Max * m_pixelsPerMeter, IM_COL32(0, 255, 0, 255), 0.0f, 0, 1.0f);
                }
                if (m_drawVelocities) {
                    auto position = origin + body->GetPosition() * m_pixelsPerMeter;
                    auto velocity = body->GetVelocity() * m_pixelsPerMeter * 0.1f;
                    drawList->AddLine(position, position + velocity, IM_COL32(0, 255, 0, 255));
                }
            }
        }
        drawList->PopClipRect();
    }
    End();
}

inline void World::ContextMenu(ImVec2 mouse) {
    using namespace ImGui;

    mouse /= m_pixelsPerMeter;
    if (MenuItem("Circle")) {
        AddBody<Circle>(1.0f, 0.5f, mouse, 1.0f);
    }
    if (MenuItem("Box")) {
        AddBody<Box>(1.0f, 0.5f, mouse, ImVec2(1.0f, 1.0f));
    }
    if (MenuItem("Poly")) {
        AddBody<Poly>(1.0f, 0.5f, mouse, std::vector<ImVec2>{
            ImVec2(-1.0f,  0.0f),
            ImVec2( 1.0f,  0.0f),
            ImVec2( 0.0f, -2.0f),
        });
    }
    Separator();
    if (MenuItem("GIMME DA BLOOST")) {
        for (float y = -5.0f; y < 5.0f; y += 1.0f) {
            for (float x = -5.0f; x < 5.0f; x += 1.0f) {
                AddBody<Circle>(1.0f, 0.5f, mouse + ImVec2(x * 0.25f, y * 0.25f), 0.25f);
            }
        }
    }
}

inline void World::DrawGrid(ImDrawList *drawList, ImRect canvasBounds, ImVec2 scrolling) {
    for (float x = ImFmod(scrolling.x, m_pixelsPerMeter); x < canvasBounds.GetWidth(); x += m_pixelsPerMeter) {
        drawList->AddLine(ImVec2(canvasBounds.Min.x + x, canvasBounds.Min.y), ImVec2(canvasBounds.Min.x + x, canvasBounds.Max.y), IM_COL32(200, 200, 200, 40));
    }
    for (float y = ImFmod(scrolling.y, m_pixelsPerMeter); y < canvasBounds.GetHeight(); y += m_pixelsPerMeter) {
        drawList->AddLine(ImVec2(canvasBounds.Min.x, canvasBounds.Min.y + y), ImVec2(canvasBounds.Max.x, canvasBounds.Min.y + y), IM_COL32(200, 200, 200, 40));
    }
}

//==-----------------------------------------------------------------------==//
// [SECTION] Entry point, more or less...
//==-----------------------------------------------------------------------==//

World world; // Yep :)

void init() {
    stm_setup();
    sg_setup({.logger = { .func = slog_func }, .environment = sglue_environment()});
    simgui_setup({.ini_filename = nullptr});
    ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    ImGui::GetStyle().AntiAliasedLines = false;

    // Easiest way to set up desired window layout
    const char *config = 
        "[Window][DockSpaceViewport_11111111]\n"
        "Size=1366,768\n"
        "Collapsed=0\n\n"
        "[Window][Viewport]\n"
        "Pos=0,0\n"
        "Size=898,768\n"
        "Collapsed=0\n"
        "DockId=0x00000001,0\n\n"
        "[Window][Config]\n"
        "Pos=900,0\n"
        "Size=466,768\n"
        "Collapsed=0\n"
        "DockId=0x00000002,0\n\n"
        "[Window][Debug##Default]\n"
        "Pos=60,60\n"
        "Size=400,400\n"
        "Collapsed=0\n\n"
        "[Window][WindowOverViewport_11111111]\n"
        "Pos=0,0\n"
        "Size=1366,768\n"
        "Collapsed=0\n\n"
        "[Docking][Data]\n"
        "DockSpace   ID=0x08BD597D Window=0x1BBC0F80 Pos=0,0 Size=1366,768 Split=X Selected=0xC450F867\n"
        "  DockNode  ID=0x00000001 Parent=0x08BD597D SizeRef=898,768 CentralNode=1 Selected=0xC450F867\n"
        "  DockNode  ID=0x00000002 Parent=0x08BD597D SizeRef=466,768 Selected=0x46F44F93\n"
        "DockSpace   ID=0x8B93E3BD Pos=0,0 Size=1366,768 CentralNode=1 Selected=0x995B0CF8\n";
    ImGui::LoadIniSettingsFromMemory(config);

    // Create a sandbox to play around in...
    world.AddBody<Box>(0.0f, 0.5f, ImVec2{0.0f, -10.0f}, ImVec2{10.0f, 0.25f});
    world.AddBody<Box>(0.0f, 0.5f, ImVec2{-10.0f, 0}, ImVec2{0.25f, 10.0f});
    world.AddBody<Box>(0.0f, 0.5f, ImVec2{10.0f, 0}, ImVec2{0.25f, 10.0f});
    world.AddBody<Box>(0.0f, 0.5f, ImVec2{0.0f, 10.0f}, ImVec2{10.0f, 0.25f});

    // ...and bunch of circles...
    for (float x = -5.0f; x < 5.0f; x += 1.0f)
        world.AddBody<Circle>(1.0f, 0.5f, ImVec2{x * 1.5f, -2.0f}, 0.75f);

    // ...boxes...
    for (float x = -5.0f; x < 5.0f; x += 1.0f)
        world.AddBody<Box>(1.0f, 0.5f, ImVec2{x * 1.25f, 0.0f},
                           ImVec2{0.5f, 0.5f});

    // ...and polygons.
    for (float x = -5.0f; x < 5.0f; x += 1.0f) {
        std::vector<ImVec2> vertices = {
            ImVec2(-1.0f, 0.0f),
            ImVec2(1.0f, 0.0f),
            ImVec2(0.0f, -2.0f),
        };
        world.AddBody<Poly>(1.0f, 0.5f, ImVec2{x * 2.0f, 3.0f}, vertices);
    }
}

void frame() {
    static uint64_t lastTime;
    static ImVec2 sandboxScrolling;
    const double deltaTime = stm_sec(stm_laptime(&lastTime));
    simgui_new_frame({
        .width = sapp_width(),
        .height = sapp_height(),
        .delta_time = deltaTime,
        .dpi_scale = sapp_dpi_scale(),
    });

    auto dockspace = ImGui::DockSpaceOverViewport();
    ImGui::SetNextWindowDockID(dockspace, ImGuiCond_Appearing);
    world.WindowViewport("Viewport", sandboxScrolling);
    world.WindowConfig("Config");
    world.Update(float(deltaTime));

    sg_begin_pass({
        .action = {
            .colors = {{
                .load_action = SG_LOADACTION_CLEAR,
                .store_action = SG_STOREACTION_STORE,
                .clear_value = {0.0f, 0.0f, 0.0f, 1.0f},
            }},
        },
        .swapchain = sglue_swapchain(),
    });
    simgui_render();
    sg_end_pass();
    sg_commit();
}

void cleanup() {
    simgui_shutdown();
    sg_shutdown();
}

void event(const sapp_event *event) {
    if (simgui_handle_event(event)) {
        return;
    }
}

sapp_desc sokol_main(int argc, char *argv[]) {
    return {
        .init_cb = init,
        .frame_cb = frame,
        .cleanup_cb = cleanup,
        .event_cb = event,
        .width = 1366,
        .height = 768,
        .sample_count = 1,
        .window_title = "rx_phys",
        .enable_clipboard = true,
    };
}
