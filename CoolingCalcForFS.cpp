#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

struct DataPoint {
    double time; // [s]
    double rpm;
    double temp_radiator_in;  // [C]
    double temp_radiator_out; // [C]
    double flow_radiator;     // [l/min]
    double temp_exchanger_in; // [C]
    double temp_exchanger_out; // [C]
    double flow_exchanger;    // [l/min]
    double temp_oil_in;       // [C]
    double temp_oil_out;      // [C]
};

struct SimulationResult {
    double time;
    double rpm;
    double radiator_cooling_power;     // [W]
    double exchanger_cooling_power;    // [W]
    double oil_heat_removed;           // [W]
    double engine_heat_generated;      // [W]
    double temp_radiator_in;           // [C]
    double temp_radiator_out;          // [C]
    double temp_exchanger_in;          // [C]
    double temp_exchanger_out;         // [C]
    double temp_oil_in;                // [C]
    double temp_oil_out;               // [C]
    bool water_overtemp;               // true if any water temp > 125C
    bool oil_overtemp;                 // true if any oil temp > 135C
    bool overheating;                  // true if engine heat > cooling capacity
};

// Physical constants
constexpr double WATER_DENSITY = 997.0;        // [kg/m³]
constexpr double WATER_SPECIFIC_HEAT = 4180.0; // [J/(kg·K)]
constexpr double OIL_DENSITY = 870.0;          // [kg/m³]
constexpr double OIL_SPECIFIC_HEAT = 2000.0;   // [J/(kg·K)]
constexpr double WATER_TEMP_LIMIT = 125.0;     // [C]
constexpr double OIL_TEMP_LIMIT = 135.0;       // [C]

// Engine power data from dynamometer (rpm, power in watts)
const std::vector<std::pair<double, double>> ENGINE_POWER_MAP = {
    {1000, 5000},
    {2000, 15000},
    {3000, 30000},
    {4000, 45000},
    {5000, 58000},
    {6000, 70000},
    {7000, 80000},
    {8000, 90000},
    {9000, 95000},
    {10000, 100000},
    {11000, 105000},
    {12000, 108000},
    {13000, 110000}
};

// Function to load log file
std::vector<DataPoint> loadLogFile(const std::string& filename) {
    std::vector<DataPoint> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // Skip header line

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        DataPoint dp;
        std::string value;

        // Parse each field from the CSV line
        std::getline(ss, value, ';'); dp.time = std::stod(value);
        std::getline(ss, value, ';'); dp.rpm = std::stod(value);
        std::getline(ss, value, ';'); dp.temp_radiator_in = std::stod(value);
        std::getline(ss, value, ';'); dp.temp_radiator_out = std::stod(value);
        std::getline(ss, value, ';'); dp.temp_exchanger_in = std::stod(value);
        std::getline(ss, value, ';'); dp.temp_exchanger_out = std::stod(value);
        std::getline(ss, value, ';'); dp.temp_oil_out = std::stod(value);
        std::getline(ss, value, ';'); dp.temp_oil_in = std::stod(value);
        std::getline(ss, value, ';'); dp.flow_radiator = std::stod(value) * 60 / 356; // Convert to l/min
        std::getline(ss, value, ';'); dp.flow_exchanger = std::stod(value) * 60 / 356; // Convert to l/min

        data.push_back(dp);
    }

    return data;
}

// Convert flow from l/min to kg/s
double lpm_to_kgps(double flow_lpm, double density) {
    // Convert from l/min to m³/s, then to kg/s
    return (flow_lpm / 60.0) * (density / 1000.0);
}

// Calculate cooling power from water flow and temperature difference
double calculate_water_cooling_power(double flow_lpm, double delta_T) {
    // P = m_dot * c * ΔT
    double flow_kgps = lpm_to_kgps(flow_lpm, WATER_DENSITY);
    return flow_kgps * WATER_SPECIFIC_HEAT * delta_T; // [W]
}



// Interpolate engine heat generation based on RPM
double interpolate_engine_heat(double rpm, const std::vector<std::pair<double, double>>& rpm_power_map) {
    // Assume 30% of engine power becomes heat
    const double HEAT_FACTOR = 0.3;

    // If RPM is below the lowest value, extrapolate from first two points
    if (rpm < rpm_power_map.front().first) {
        double rpm1 = rpm_power_map[0].first;
        double rpm2 = rpm_power_map[1].first;
        double power1 = rpm_power_map[0].second;
        double power2 = rpm_power_map[1].second;
        double power = power1 + (rpm - rpm1) * (power2 - power1) / (rpm2 - rpm1);
        return std::max(0.0, power * HEAT_FACTOR);
    }

    // If RPM is above the highest value, use the highest value
    if (rpm >= rpm_power_map.back().first) {
        return rpm_power_map.back().second * HEAT_FACTOR;
    }

    // Otherwise, interpolate between two closest points
    for (size_t i = 0; i < rpm_power_map.size() - 1; ++i) {
        if (rpm >= rpm_power_map[i].first && rpm < rpm_power_map[i + 1].first) {
            double rpm1 = rpm_power_map[i].first;
            double rpm2 = rpm_power_map[i + 1].first;
            double power1 = rpm_power_map[i].second;
            double power2 = rpm_power_map[i + 1].second;
            double power = power1 + (rpm - rpm1) * (power2 - power1) / (rpm2 - rpm1);
            return power * HEAT_FACTOR;
        }
    }

    // This should not happen, but return a default value just in case
    return rpm_power_map.back().second * HEAT_FACTOR;
}

// Calculate radiator effectiveness based on ambient temperature
double calculate_radiator_effectiveness(double ambient_temp, double original_ambient_temp = 25.0) {
    // As ambient temperature approaches the coolant temperature, effectiveness decreases
    // This is a simplified model based on heat transfer theory

    // Base effectiveness at original ambient temperature
    const double BASE_EFFECTIVENESS = 1.0;

    // Maximum effectiveness reduction at very high ambient temperatures
    const double MAX_EFFECTIVENESS_REDUCTION = 0.5;

    // Temperature difference at which effectiveness starts to significantly decline
    const double EFFECTIVENESS_DECLINE_THRESHOLD = 20.0;

    // Calculate temperature increase from baseline
    double temp_increase = ambient_temp - original_ambient_temp;

    // If temperature is below or equal to original, no reduction in effectiveness
    if (temp_increase <= 0.0) {
        return BASE_EFFECTIVENESS;
    }

    // Calculate effectiveness reduction factor
    double reduction_factor = std::min(1.0, temp_increase / EFFECTIVENESS_DECLINE_THRESHOLD);

    // Calculate new effectiveness
    double effectiveness = BASE_EFFECTIVENESS - (MAX_EFFECTIVENESS_REDUCTION * reduction_factor);

    return std::max(0.2, effectiveness); // Never allow effectiveness to go below 20%
}

// Calculate realistic temperature changes based on ambient temperature increase
DataPoint adjust_for_ambient_temp(const DataPoint& dp, double ambient_temp, double original_ambient_temp = 25.0) {
    DataPoint adjusted = dp;

    // Calculate radiator effectiveness based on new ambient temperature
    double radiator_effectiveness = calculate_radiator_effectiveness(ambient_temp, original_ambient_temp);

    // Calculate temperature difference between water and ambient in normal conditions
    double normal_radiator_water_ambient_diff = dp.temp_radiator_out - original_ambient_temp;

    // Calculate new radiator outlet temperature based on reduced effectiveness
    // The closer the outlet temperature gets to ambient, the harder it is to cool further
    double new_radiator_outlet_diff = normal_radiator_water_ambient_diff * radiator_effectiveness;
    adjusted.temp_radiator_out = ambient_temp + new_radiator_outlet_diff;

    // Calculate normal temperature drop across radiator
    double normal_radiator_drop = dp.temp_radiator_in - dp.temp_radiator_out;

    // Apply same drop to calculate new inlet temperature (with some reduction due to higher overall temps)
    adjusted.temp_radiator_in = adjusted.temp_radiator_out + (normal_radiator_drop * radiator_effectiveness);

    // Calculate normal temperature difference between radiator inlet and exchanger outlet
    double normal_exchanger_out_diff = dp.temp_exchanger_out - dp.temp_radiator_in;

    // Apply same difference to calculate new exchanger outlet temperature
    adjusted.temp_exchanger_out = adjusted.temp_radiator_in + normal_exchanger_out_diff;

    // Calculate normal temperature drop across exchanger water side
    double normal_exchanger_water_drop = dp.temp_exchanger_in - dp.temp_exchanger_out;

    // Apply same drop to calculate new exchanger water inlet temperature
    adjusted.temp_exchanger_in = adjusted.temp_exchanger_out + normal_exchanger_water_drop;

    // Calculate effectiveness of oil cooler based on higher water temperatures
    double oil_cooler_effectiveness = radiator_effectiveness * 0.9; // Slightly worse than radiator

    // Calculate normal temperature difference between oil in and water in at exchanger
    double normal_oil_water_diff = dp.temp_oil_in - dp.temp_exchanger_in;

    // Calculate new oil inlet temperature based on new water inlet temperature and reduced ability to cool
    adjusted.temp_oil_in = adjusted.temp_exchanger_in + (normal_oil_water_diff / oil_cooler_effectiveness);

    // Calculate normal temperature drop across oil side of exchanger
    double normal_oil_drop = dp.temp_oil_in - dp.temp_oil_out;

    // Apply same drop to calculate new oil outlet temperature (with reduction due to less effective cooling)
    adjusted.temp_oil_out = adjusted.temp_oil_in - (normal_oil_drop * oil_cooler_effectiveness);

    return adjusted;
}

// Run simulation with given parameters
void simulate(const std::vector<DataPoint>& data,
    std::vector<SimulationResult>& results,
    double ambient_temp = 25.0,
    double original_ambient_temp = 25.0,
    double extra_time = 0.0) {

    // Clear previous results
    results.clear();

    // Process each data point
    for (size_t i = 0; i < data.size(); ++i) {
        // Calculate time step (dt)
        double dt = (i > 0) ? (data[i].time - data[i - 1].time) : 0.04; // Assume 25Hz = 0.04s if first point

        // Adjust temperatures for ambient increase
        DataPoint adjusted = adjust_for_ambient_temp(data[i], ambient_temp, original_ambient_temp);

        // Calculate temperature differences
        double delta_T_radiator = adjusted.temp_radiator_in - adjusted.temp_radiator_out;
        double delta_T_exchanger_water = adjusted.temp_exchanger_out - adjusted.temp_exchanger_in;
        double delta_T_oil = adjusted.temp_oil_in - adjusted.temp_oil_out;

        // Calculate cooling powers
        // Note: Multiply radiator flow by 2 because there are two radiators
        double radiator_cooling_power = calculate_water_cooling_power(2 * adjusted.flow_radiator, delta_T_radiator);
        double exchanger_cooling_power = calculate_water_cooling_power(adjusted.flow_exchanger, delta_T_exchanger_water);

        // Calculate oil heat removed (estimated from water side of exchanger)
        double oil_heat_removed = exchanger_cooling_power;

        // Estimate engine heat generation from RPM
        double engine_heat = interpolate_engine_heat(adjusted.rpm, ENGINE_POWER_MAP);

        // Check temperature thresholds
        bool water_overtemp = (adjusted.temp_radiator_in > WATER_TEMP_LIMIT) ||
            (adjusted.temp_radiator_out > WATER_TEMP_LIMIT) ||
            (adjusted.temp_exchanger_in > WATER_TEMP_LIMIT) ||
            (adjusted.temp_exchanger_out > WATER_TEMP_LIMIT);

        bool oil_overtemp = (adjusted.temp_oil_in > OIL_TEMP_LIMIT) ||
            (adjusted.temp_oil_out > OIL_TEMP_LIMIT);

        // Check if cooling is sufficient
        bool overheating = engine_heat > oil_heat_removed;

        // Store results
        SimulationResult res;
        res.time = adjusted.time;
        res.rpm = adjusted.rpm;
        res.radiator_cooling_power = radiator_cooling_power;
        res.exchanger_cooling_power = exchanger_cooling_power;
        res.oil_heat_removed = oil_heat_removed;
        res.engine_heat_generated = engine_heat;
        res.temp_radiator_in = adjusted.temp_radiator_in;
        res.temp_radiator_out = adjusted.temp_radiator_out;
        res.temp_exchanger_in = adjusted.temp_exchanger_in;
        res.temp_exchanger_out = adjusted.temp_exchanger_out;
        res.temp_oil_in = adjusted.temp_oil_in;
        res.temp_oil_out = adjusted.temp_oil_out;
        res.water_overtemp = water_overtemp;
        res.oil_overtemp = oil_overtemp;
        res.overheating = overheating;

        results.push_back(res);
    }

    // Handle extended simulation time if needed
    if (extra_time > 0.0 && !data.empty()) {
        double end_time = data.back().time;
        double target_end_time = end_time + extra_time;
        double time_step = 0.04; // 25Hz

        // Use the last data point as a base for extended simulation
        DataPoint last_point = data.back();
        SimulationResult last_result = results.back();

        // Simulation parameters for continuous operation
        double current_engine_heat = last_result.engine_heat_generated;
        double current_exchanger_power = last_result.exchanger_cooling_power;
        double thermal_imbalance = current_engine_heat - current_exchanger_power;

        // Temperature rise rates (°C per second) if there's an imbalance
        // These are simplified models - in real life would need thermal mass calculations
        const double WATER_TEMP_RISE_RATE = 0.05; // °C per second per kW of imbalance
        const double OIL_TEMP_RISE_RATE = 0.08; // °C per second per kW of imbalance

        // Add extended time data points
        for (double t = end_time + time_step; t <= target_end_time; t += time_step) {
            SimulationResult res = last_result;
            res.time = t;

            // If there's a cooling deficit, temperatures will rise
            if (thermal_imbalance > 0) {
                double time_elapsed = t - end_time;
                double temp_increase = (thermal_imbalance / 1000.0) * WATER_TEMP_RISE_RATE * time_elapsed;
                double oil_temp_increase = (thermal_imbalance / 1000.0) * OIL_TEMP_RISE_RATE * time_elapsed;

                // Apply temperature increases
                res.temp_radiator_in += temp_increase;
                res.temp_radiator_out += temp_increase;
                res.temp_exchanger_in += temp_increase;
                res.temp_exchanger_out += temp_increase;
                res.temp_oil_in += oil_temp_increase;
                res.temp_oil_out += oil_temp_increase;

                // Recalculate cooling powers with new temperatures
                double delta_T_radiator = res.temp_radiator_in - res.temp_radiator_out;
                double delta_T_exchanger = res.temp_exchanger_out - res.temp_exchanger_in;

                // Use the same flow rates as the last data point
                res.radiator_cooling_power = calculate_water_cooling_power(2 * last_point.flow_radiator, delta_T_radiator);
                res.exchanger_cooling_power = calculate_water_cooling_power(last_point.flow_exchanger, delta_T_exchanger);

                // As temperatures rise, cooling effectiveness increases (bigger ΔT with ambient)
                // This creates a balancing effect
                double cooling_increase_factor = 1.0 + (temp_increase / 100.0); // Simple model
                res.radiator_cooling_power *= cooling_increase_factor;
                res.exchanger_cooling_power *= cooling_increase_factor;
                res.oil_heat_removed = res.exchanger_cooling_power;
            }

            // Update overheating flags
            res.water_overtemp = (res.temp_radiator_in > WATER_TEMP_LIMIT) ||
                (res.temp_radiator_out > WATER_TEMP_LIMIT) ||
                (res.temp_exchanger_in > WATER_TEMP_LIMIT) ||
                (res.temp_exchanger_out > WATER_TEMP_LIMIT);

            res.oil_overtemp = (res.temp_oil_in > OIL_TEMP_LIMIT) ||
                (res.temp_oil_out > OIL_TEMP_LIMIT);

            res.overheating = res.engine_heat_generated > res.oil_heat_removed;

            results.push_back(res);
        }
    }
}

// Save simulation results to CSV
void save_simulation(const std::string& filename, const std::vector<SimulationResult>& results) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to create output file: " << filename << std::endl;
        return;
    }

    // Write CSV header
    file << "time;rpm;radiator_cooling_power;exchanger_cooling_power;oil_heat_removed;"
        << "engine_heat_generated;temp_radiator_in;temp_radiator_out;temp_exchanger_in;"
        << "temp_exchanger_out;temp_oil_in;temp_oil_out;water_overtemp;oil_overtemp;overheating\n";

    // Write data rows
    for (const auto& res : results) {
        file << std::fixed << std::setprecision(3)
            << res.time << ";"
            << res.rpm << ";"
            << res.radiator_cooling_power << ";"
            << res.exchanger_cooling_power << ";"
            << res.oil_heat_removed << ";"
            << res.engine_heat_generated << ";"
            << res.temp_radiator_in << ";"
            << res.temp_radiator_out << ";"
            << res.temp_exchanger_in << ";"
            << res.temp_exchanger_out << ";"
            << res.temp_oil_in << ";"
            << res.temp_oil_out << ";"
            << (res.water_overtemp ? 1 : 0) << ";"
            << (res.oil_overtemp ? 1 : 0) << ";"
            << (res.overheating ? 1 : 0) << "\n";
    }

    file.close();
}

// Calculate statistical summary of simulation results
void analyze_simulation(const std::vector<SimulationResult>& results, const std::string& scenario_name) {
    if (results.empty()) {
        std::cout << "No results to analyze for scenario: " << scenario_name << std::endl;
        return;
    }

    // Calculate statistics
    double avg_radiator_power = 0.0;
    double avg_exchanger_power = 0.0;
    double avg_engine_heat = 0.0;
    double avg_water_in_temp = 0.0;
    double avg_oil_in_temp = 0.0;
    double max_water_temp = 0.0;
    double max_oil_temp = 0.0;
    int overheating_count = 0;
    int water_overtemp_count = 0;
    int oil_overtemp_count = 0;

    for (const auto& res : results) {
        avg_radiator_power += res.radiator_cooling_power;
        avg_exchanger_power += res.exchanger_cooling_power;
        avg_engine_heat += res.engine_heat_generated;
        avg_water_in_temp += res.temp_radiator_in;
        avg_oil_in_temp += res.temp_oil_in;

        max_water_temp = std::max({ max_water_temp, res.temp_radiator_in, res.temp_radiator_out,
                                  res.temp_exchanger_in, res.temp_exchanger_out });
        max_oil_temp = std::max({ max_oil_temp, res.temp_oil_in, res.temp_oil_out });

        if (res.overheating) overheating_count++;
        if (res.water_overtemp) water_overtemp_count++;
        if (res.oil_overtemp) oil_overtemp_count++;
    }

    // Calculate averages
    avg_radiator_power /= results.size();
    avg_exchanger_power /= results.size();
    avg_engine_heat /= results.size();
    avg_water_in_temp /= results.size();
    avg_oil_in_temp /= results.size();

    double overheating_percent = 100.0 * overheating_count / results.size();
    double water_overtemp_percent = 100.0 * water_overtemp_count / results.size();
    double oil_overtemp_percent = 100.0 * oil_overtemp_count / results.size();

    // Print results
    std::cout << "=== Analysis for " << scenario_name << " ===" << std::endl;
    std::cout << "Average radiator cooling power: " << avg_radiator_power << " W" << std::endl;
    std::cout << "Average exchanger cooling power: " << avg_exchanger_power << " W" << std::endl;
    std::cout << "Average engine heat generation: " << avg_engine_heat << " W" << std::endl;
    std::cout << "Average water inlet temperature: " << avg_water_in_temp << " °C" << std::endl;
    std::cout << "Average oil inlet temperature: " << avg_oil_in_temp << " °C" << std::endl;
    std::cout << "Maximum water temperature: " << max_water_temp << " °C" << std::endl;
    std::cout << "Maximum oil temperature: " << max_oil_temp << " °C" << std::endl;
    std::cout << "Overheating occurred in " << overheating_count << " samples ("
        << overheating_percent << "%)" << std::endl;
    std::cout << "Water over " << WATER_TEMP_LIMIT << "°C in " << water_overtemp_count << " samples ("
        << water_overtemp_percent << "%)" << std::endl;
    std::cout << "Oil over " << OIL_TEMP_LIMIT << "°C in " << oil_overtemp_count << " samples ("
        << oil_overtemp_percent << "%)" << std::endl;
    std::cout << std::endl;
}

int main() {
    // Configuration
    std::string input_file = "przejazdMarekZiecLogAdu.txt";
    double original_ambient_temp = 25.0; // Assumed ambient temperature during original data collection

    std::cout << "Loading data from " << input_file << "..." << std::endl;
    auto data = loadLogFile(input_file);

    if (data.empty()) {
        std::cerr << "No data loaded. Exiting." << std::endl;
        return 1;
    }

    std::cout << "Loaded " << data.size() << " data points." << std::endl;

    // Run simulations with different scenarios

    // Scenario 1: Current conditions (baseline)
    std::vector<SimulationResult> normal_results;
    std::cout << "Running simulation with current conditions..." << std::endl;
    simulate(data, normal_results, original_ambient_temp, original_ambient_temp, 0.0);
    save_simulation("simulation_current.csv", normal_results);
    analyze_simulation(normal_results, "Current Conditions");

    // Scenario 2: Increased ambient temperature (+10°C)
    std::vector<SimulationResult> hot_weather_results;
    std::cout << "Running simulation with +10°C ambient temperature..." << std::endl;
    simulate(data, hot_weather_results, original_ambient_temp + 10.0, original_ambient_temp, 0.0);
    save_simulation("simulation_hot.csv", hot_weather_results);
    analyze_simulation(hot_weather_results, "Hot Weather (+10°C)");

    // Scenario 3: Increased ambient temperature (+10°C) with extended run time (+600s)
    std::vector<SimulationResult> extended_hot_results;
    std::cout << "Running simulation with +10°C ambient temperature and extended time..." << std::endl;
    simulate(data, extended_hot_results, original_ambient_temp + 10.0, original_ambient_temp, 600.0);
    save_simulation("simulation_hot_extended.csv", extended_hot_results);
    analyze_simulation(extended_hot_results, "Extended Hot Weather (+10°C, +600s)");

    // Scenario 4: Extreme conditions (+20°C)
    std::vector<SimulationResult> extreme_hot_results;
    std::cout << "Running simulation with +20°C ambient temperature..." << std::endl;
    simulate(data, extreme_hot_results, original_ambient_temp + 20.0, original_ambient_temp, 0.0);
    save_simulation("simulation_extreme_hot.csv", extreme_hot_results);
    analyze_simulation(extreme_hot_results, "Extreme Heat (+20°C)");

    std::cout << "All simulations completed." << std::endl;
    return 0;
}