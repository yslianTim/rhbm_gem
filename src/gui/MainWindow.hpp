#pragma once

#include <QMainWindow>
#include <QString>

#include <future>

#include <rhbm_gem/core/command/CommandApi.hpp>

class QCheckBox;
class QComboBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QSpinBox;
class QDoubleSpinBox;
class QStackedWidget;
class QFormLayout;
class QTimer;

class MainWindow final : public QMainWindow
{
public:
    explicit MainWindow(QWidget * parent = nullptr);
    ~MainWindow() override = default;

private:
    struct CommonControls
    {
        QSpinBox * jobs{ nullptr };
        QSpinBox * verbose{ nullptr };
        QLineEdit * database_path{ nullptr };
        QPushButton * database_browse{ nullptr };
        QLineEdit * folder_path{ nullptr };
        QPushButton * folder_browse{ nullptr };
    };

    struct MapSimulationControls
    {
        CommonControls common{};
        QLineEdit * model_path{ nullptr };
        QPushButton * model_browse{ nullptr };
        QLineEdit * map_name{ nullptr };
        QComboBox * potential_model{ nullptr };
        QComboBox * partial_charge{ nullptr };
        QDoubleSpinBox * cutoff_distance{ nullptr };
        QDoubleSpinBox * grid_spacing{ nullptr };
        QLineEdit * blurring_widths{ nullptr };
    };

    struct PotentialAnalysisControls
    {
        CommonControls common{};
        QLineEdit * model_path{ nullptr };
        QPushButton * model_browse{ nullptr };
        QLineEdit * map_path{ nullptr };
        QPushButton * map_browse{ nullptr };
        QLineEdit * saved_key_tag{ nullptr };
        QLineEdit * training_report_dir{ nullptr };
        QPushButton * training_report_dir_browse{ nullptr };
        QCheckBox * simulation_flag{ nullptr };
        QDoubleSpinBox * simulated_map_resolution{ nullptr };
        QCheckBox * training_alpha_flag{ nullptr };
        QCheckBox * asymmetry_flag{ nullptr };
        QSpinBox * sampling_size{ nullptr };
        QDoubleSpinBox * sampling_range_min{ nullptr };
        QDoubleSpinBox * sampling_range_max{ nullptr };
        QDoubleSpinBox * sampling_height{ nullptr };
        QDoubleSpinBox * fit_range_min{ nullptr };
        QDoubleSpinBox * fit_range_max{ nullptr };
        QDoubleSpinBox * alpha_r{ nullptr };
        QDoubleSpinBox * alpha_g{ nullptr };
    };

    struct ResultDumpControls
    {
        CommonControls common{};
        QComboBox * printer_choice{ nullptr };
        QLineEdit * model_key_list{ nullptr };
        QLineEdit * map_path{ nullptr };
        QPushButton * map_browse{ nullptr };
    };

    enum class CommandPage
    {
        MapSimulation = 0,
        PotentialAnalysis = 1,
        ResultDump = 2
    };

    QListWidget * m_command_list{ nullptr };
    QStackedWidget * m_command_stack{ nullptr };
    QPushButton * m_run_button{ nullptr };
    QLabel * m_status_label{ nullptr };
    QPlainTextEdit * m_result_view{ nullptr };
    QTimer * m_execution_poll_timer{ nullptr };
    bool m_execution_running{ false };
    std::future<rhbm_gem::ExecutionReport> m_execution_future{};
    QString m_active_command_name{};

    MapSimulationControls m_map_simulation{};
    PotentialAnalysisControls m_potential_analysis{};
    ResultDumpControls m_result_dump{};

    void BuildUi();
    QWidget * BuildMapSimulationPage();
    QWidget * BuildPotentialAnalysisPage();
    QWidget * BuildResultDumpPage();
    void AddCommonControls(
        QFormLayout * layout,
        CommonControls & controls,
        bool include_database);
    void ConnectUi();

    void BrowseForInputFile(QLineEdit * target, const QString & caption);
    void BrowseForOutputFile(QLineEdit * target, const QString & caption);
    void BrowseForDirectory(QLineEdit * target, const QString & caption);

    void StartExecution();
    void PollExecutionState();
    void OnExecutionFinished(const rhbm_gem::ExecutionReport & result);

    rhbm_gem::MapSimulationRequest BuildMapSimulationRequest() const;
    rhbm_gem::PotentialAnalysisRequest BuildPotentialAnalysisRequest() const;
    rhbm_gem::ResultDumpRequest BuildResultDumpRequest() const;

    static QString FormatExecutionSummary(const rhbm_gem::ExecutionReport & result);
    static QString FormatValidationIssues(
        const std::vector<rhbm_gem::ValidationIssue> & issues);
    static QString ValidationPhaseLabel(rhbm_gem::ValidationPhase phase);
    static QString ValidationLevelLabel(LogLevel level);
};
