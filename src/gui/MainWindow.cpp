#include "MainWindow.hpp"

#include "internal/CommandCatalog.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QDateTime>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSplitter>
#include <QSpinBox>
#include <QStackedWidget>
#include <QTimer>
#include <QVBoxLayout>
#include <QWidget>

#include <rhbm_gem/core/command/OptionEnumTraits.hpp>

#include <algorithm>
#include <chrono>
#include <stdexcept>

namespace {

const rhbm_gem::CommandDescriptor & RequireCommandDescriptor(rhbm_gem::CommandId command_id)
{
    const auto & catalog{ rhbm_gem::CommandCatalog() };
    const auto iter{
        std::find_if(
            catalog.begin(),
            catalog.end(),
            [command_id](const rhbm_gem::CommandDescriptor & descriptor)
            {
                return descriptor.id == command_id;
            })
    };
    if (iter == catalog.end())
    {
        throw std::runtime_error("GUI command descriptor is missing from CommandCatalog().");
    }
    return *iter;
}

QString CommandName(rhbm_gem::CommandId command_id)
{
    const auto & descriptor{ RequireCommandDescriptor(command_id) };
    return QString::fromUtf8(descriptor.name.data(), static_cast<int>(descriptor.name.size()));
}

bool CommandUsesDatabase(rhbm_gem::CommandId command_id)
{
    return rhbm_gem::HasCommonOption(
        rhbm_gem::CommonOptionsForCommand(RequireCommandDescriptor(command_id)),
        rhbm_gem::CommonOption::Database);
}

QWidget * BuildPathSelector(QLineEdit *& line_edit, QPushButton *& browse_button)
{
    auto * row{ new QWidget() };
    auto * row_layout{ new QHBoxLayout(row) };
    row_layout->setContentsMargins(0, 0, 0, 0);
    row_layout->setSpacing(8);

    line_edit = new QLineEdit(row);
    browse_button = new QPushButton("Browse...", row);

    row_layout->addWidget(line_edit, 1);
    row_layout->addWidget(browse_button);
    return row;
}

template <typename EnumType>
void PopulateEnumCombo(QComboBox * combo)
{
    const auto binding_entries{ rhbm_gem::GetEnumBindingEntries<EnumType>() };
    for (const auto & entry : binding_entries)
    {
        const QString label{
            QString::fromUtf8(entry.token.data(), static_cast<int>(entry.token.size()))
        };
        combo->addItem(label, static_cast<int>(entry.value));
    }
}

template <typename EnumType>
EnumType CurrentEnumValue(const QComboBox * combo)
{
    return static_cast<EnumType>(combo->currentData().toInt());
}

void SetCurrentEnumValue(QComboBox * combo, int value)
{
    const int index{ combo->findData(value) };
    if (index >= 0)
    {
        combo->setCurrentIndex(index);
    }
}

std::filesystem::path ReadPath(const QLineEdit * line_edit)
{
    const auto path_utf8{ line_edit->text().trimmed().toUtf8() };
    return std::filesystem::path(
        std::string(path_utf8.constData(), static_cast<std::size_t>(path_utf8.size())));
}

std::string ReadString(const QLineEdit * line_edit)
{
    const auto text_utf8{ line_edit->text().trimmed().toUtf8() };
    return std::string(
        text_utf8.constData(), static_cast<std::size_t>(text_utf8.size()));
}

} // namespace

MainWindow::MainWindow(QWidget * parent) :
    QMainWindow(parent)
{
    BuildUi();
    ConnectUi();
}

void MainWindow::InitializeGuiCommands()
{
    if (!m_gui_commands.empty())
    {
        return;
    }

    m_gui_commands = {
        GuiCommandRegistration{
            rhbm_gem::CommandId::MapSimulation,
            [this]() { return BuildMapSimulationPage(); },
            [this]() { return rhbm_gem::RunMapSimulation(BuildMapSimulationRequest()); }
        },
        GuiCommandRegistration{
            rhbm_gem::CommandId::PotentialAnalysis,
            [this]() { return BuildPotentialAnalysisPage(); },
            [this]() { return rhbm_gem::RunPotentialAnalysis(BuildPotentialAnalysisRequest()); }
        },
        GuiCommandRegistration{
            rhbm_gem::CommandId::ResultDump,
            [this]() { return BuildResultDumpPage(); },
            [this]() { return rhbm_gem::RunResultDump(BuildResultDumpRequest()); }
        },
    };
}

void MainWindow::BuildUi()
{
    InitializeGuiCommands();

    setWindowTitle("RHBM-GEM GUI");
    resize(1180, 820);

    auto * central{ new QWidget(this) };
    auto * root_layout{ new QVBoxLayout(central) };
    root_layout->setContentsMargins(12, 12, 12, 12);
    root_layout->setSpacing(10);

    auto * top_splitter{ new QSplitter(Qt::Horizontal, central) };
    m_command_list = new QListWidget(top_splitter);
    m_command_list->setMinimumWidth(220);

    m_command_stack = new QStackedWidget(top_splitter);
    for (const auto & command : m_gui_commands)
    {
        m_command_list->addItem(CommandName(command.id));
        m_command_stack->addWidget(command.build_page());
    }

    top_splitter->addWidget(m_command_list);
    top_splitter->addWidget(m_command_stack);
    top_splitter->setStretchFactor(0, 0);
    top_splitter->setStretchFactor(1, 1);

    auto * action_layout{ new QHBoxLayout() };
    m_run_button = new QPushButton("Run", central);
    m_status_label = new QLabel("Idle", central);
    action_layout->addWidget(m_run_button);
    action_layout->addWidget(m_status_label, 1);

    m_result_view = new QPlainTextEdit(central);
    m_result_view->setReadOnly(true);
    m_result_view->setMinimumHeight(220);
    m_result_view->setPlaceholderText(
        "Execution results and validation diagnostics will appear here.");

    root_layout->addWidget(top_splitter, 1);
    root_layout->addLayout(action_layout);
    root_layout->addWidget(m_result_view);

    setCentralWidget(central);
    m_command_list->setCurrentRow(0);
}

QWidget * MainWindow::BuildMapSimulationPage()
{
    auto * page{ new QWidget(this) };
    auto * layout{ new QFormLayout(page) };
    layout->setLabelAlignment(Qt::AlignLeft);
    layout->setFormAlignment(Qt::AlignTop);
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    AddCommonControls(
        layout,
        m_map_simulation.common,
        CommandUsesDatabase(rhbm_gem::CommandId::MapSimulation));

    layout->addRow("Model File (--model)", BuildPathSelector(
        m_map_simulation.model_path, m_map_simulation.model_browse));

    m_map_simulation.map_name = new QLineEdit("sim_map", page);
    layout->addRow("Output Map Name (--name)", m_map_simulation.map_name);

    m_map_simulation.potential_model = new QComboBox(page);
    PopulateEnumCombo<rhbm_gem::PotentialModel>(m_map_simulation.potential_model);
    SetCurrentEnumValue(
        m_map_simulation.potential_model,
        static_cast<int>(rhbm_gem::PotentialModel::FIVE_GAUS_CHARGE));
    layout->addRow("Potential Model (--potential-model)", m_map_simulation.potential_model);

    m_map_simulation.partial_charge = new QComboBox(page);
    PopulateEnumCombo<rhbm_gem::PartialCharge>(m_map_simulation.partial_charge);
    SetCurrentEnumValue(
        m_map_simulation.partial_charge,
        static_cast<int>(rhbm_gem::PartialCharge::PARTIAL));
    layout->addRow("Partial Charge (--charge)", m_map_simulation.partial_charge);

    m_map_simulation.cutoff_distance = new QDoubleSpinBox(page);
    m_map_simulation.cutoff_distance->setRange(0.0001, 1000.0);
    m_map_simulation.cutoff_distance->setDecimals(4);
    m_map_simulation.cutoff_distance->setValue(5.0);
    layout->addRow("Cutoff Distance (--cut-off)", m_map_simulation.cutoff_distance);

    m_map_simulation.grid_spacing = new QDoubleSpinBox(page);
    m_map_simulation.grid_spacing->setRange(0.0001, 1000.0);
    m_map_simulation.grid_spacing->setDecimals(4);
    m_map_simulation.grid_spacing->setValue(0.5);
    layout->addRow("Grid Spacing (--grid-spacing)", m_map_simulation.grid_spacing);

    m_map_simulation.blurring_widths = new QLineEdit("1.50", page);
    layout->addRow("Blurring Widths (--blurring-width)", m_map_simulation.blurring_widths);

    return page;
}

QWidget * MainWindow::BuildPotentialAnalysisPage()
{
    auto * page{ new QWidget(this) };
    auto * layout{ new QFormLayout(page) };
    layout->setLabelAlignment(Qt::AlignLeft);
    layout->setFormAlignment(Qt::AlignTop);
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    AddCommonControls(
        layout,
        m_potential_analysis.common,
        CommandUsesDatabase(rhbm_gem::CommandId::PotentialAnalysis));

    layout->addRow("Model File (--model)", BuildPathSelector(
        m_potential_analysis.model_path, m_potential_analysis.model_browse));
    layout->addRow("Map File (--map)", BuildPathSelector(
        m_potential_analysis.map_path, m_potential_analysis.map_browse));

    m_potential_analysis.saved_key_tag = new QLineEdit("model", page);
    layout->addRow("Saved Key (--save-key)", m_potential_analysis.saved_key_tag);
    layout->addRow(
        "Training Report Dir (--training-report-dir)",
        BuildPathSelector(
            m_potential_analysis.training_report_dir,
            m_potential_analysis.training_report_dir_browse));

    m_potential_analysis.simulation_flag = new QCheckBox(page);
    layout->addRow("Simulation Mode (--simulation)", m_potential_analysis.simulation_flag);

    m_potential_analysis.simulated_map_resolution = new QDoubleSpinBox(page);
    m_potential_analysis.simulated_map_resolution->setRange(0.0, 1000.0);
    m_potential_analysis.simulated_map_resolution->setDecimals(4);
    m_potential_analysis.simulated_map_resolution->setValue(0.0);
    m_potential_analysis.simulated_map_resolution->setEnabled(false);
    layout->addRow(
        "Simulated Map Resolution (--sim-resolution)",
        m_potential_analysis.simulated_map_resolution);

    m_potential_analysis.training_alpha_flag = new QCheckBox(page);
    layout->addRow("Training Alpha (--training-alpha)", m_potential_analysis.training_alpha_flag);

    m_potential_analysis.asymmetry_flag = new QCheckBox(page);
    layout->addRow("Asymmetry (--asymmetry)", m_potential_analysis.asymmetry_flag);

    m_potential_analysis.sampling_size = new QSpinBox(page);
    m_potential_analysis.sampling_size->setRange(1, 2000000);
    m_potential_analysis.sampling_size->setValue(1500);
    layout->addRow("Sampling Size (--sampling)", m_potential_analysis.sampling_size);

    m_potential_analysis.sampling_range_min = new QDoubleSpinBox(page);
    m_potential_analysis.sampling_range_min->setRange(0.0, 1000.0);
    m_potential_analysis.sampling_range_min->setDecimals(4);
    m_potential_analysis.sampling_range_min->setValue(0.0);
    layout->addRow("Sampling Min (--sampling-min)", m_potential_analysis.sampling_range_min);

    m_potential_analysis.sampling_range_max = new QDoubleSpinBox(page);
    m_potential_analysis.sampling_range_max->setRange(0.0, 1000.0);
    m_potential_analysis.sampling_range_max->setDecimals(4);
    m_potential_analysis.sampling_range_max->setValue(1.5);
    layout->addRow("Sampling Max (--sampling-max)", m_potential_analysis.sampling_range_max);

    m_potential_analysis.sampling_height = new QDoubleSpinBox(page);
    m_potential_analysis.sampling_height->setRange(0.0001, 1000.0);
    m_potential_analysis.sampling_height->setDecimals(4);
    m_potential_analysis.sampling_height->setValue(0.1);
    layout->addRow("Sampling Height (--sampling-height)", m_potential_analysis.sampling_height);

    m_potential_analysis.fit_range_min = new QDoubleSpinBox(page);
    m_potential_analysis.fit_range_min->setRange(0.0, 1000.0);
    m_potential_analysis.fit_range_min->setDecimals(4);
    m_potential_analysis.fit_range_min->setValue(0.0);
    layout->addRow("Fit Min (--fit-min)", m_potential_analysis.fit_range_min);

    m_potential_analysis.fit_range_max = new QDoubleSpinBox(page);
    m_potential_analysis.fit_range_max->setRange(0.0, 1000.0);
    m_potential_analysis.fit_range_max->setDecimals(4);
    m_potential_analysis.fit_range_max->setValue(1.0);
    layout->addRow("Fit Max (--fit-max)", m_potential_analysis.fit_range_max);

    m_potential_analysis.alpha_r = new QDoubleSpinBox(page);
    m_potential_analysis.alpha_r->setRange(0.0001, 1000.0);
    m_potential_analysis.alpha_r->setDecimals(6);
    m_potential_analysis.alpha_r->setValue(0.1);
    layout->addRow("Alpha R (--alpha-r)", m_potential_analysis.alpha_r);

    m_potential_analysis.alpha_g = new QDoubleSpinBox(page);
    m_potential_analysis.alpha_g->setRange(0.0001, 1000.0);
    m_potential_analysis.alpha_g->setDecimals(6);
    m_potential_analysis.alpha_g->setValue(0.2);
    layout->addRow("Alpha G (--alpha-g)", m_potential_analysis.alpha_g);

    return page;
}

QWidget * MainWindow::BuildResultDumpPage()
{
    auto * page{ new QWidget(this) };
    auto * layout{ new QFormLayout(page) };
    layout->setLabelAlignment(Qt::AlignLeft);
    layout->setFormAlignment(Qt::AlignTop);
    layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);

    AddCommonControls(
        layout,
        m_result_dump.common,
        CommandUsesDatabase(rhbm_gem::CommandId::ResultDump));

    m_result_dump.printer_choice = new QComboBox(page);
    PopulateEnumCombo<rhbm_gem::PrinterType>(m_result_dump.printer_choice);
    SetCurrentEnumValue(
        m_result_dump.printer_choice,
        static_cast<int>(rhbm_gem::PrinterType::GAUS_ESTIMATES));
    layout->addRow("Printer (--printer)", m_result_dump.printer_choice);

    m_result_dump.model_key_list = new QLineEdit(page);
    layout->addRow("Model Keys (--model-keylist)", m_result_dump.model_key_list);

    layout->addRow("Map File (--map)", BuildPathSelector(
        m_result_dump.map_path, m_result_dump.map_browse));

    const auto update_map_controls = [this]()
    {
        const bool map_required{
            CurrentEnumValue<rhbm_gem::PrinterType>(m_result_dump.printer_choice)
            == rhbm_gem::PrinterType::MAP_VALUE
        };
        m_result_dump.map_path->setEnabled(map_required);
        m_result_dump.map_browse->setEnabled(map_required);
    };

    connect(
        m_result_dump.printer_choice,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        this,
        [update_map_controls](int) { update_map_controls(); });
    update_map_controls();

    return page;
}

void MainWindow::AddCommonControls(
    QFormLayout * layout,
    CommonControls & controls,
    bool include_database)
{
    controls.jobs = new QSpinBox(this);
    controls.jobs->setRange(1, 512);
    controls.jobs->setValue(1);
    layout->addRow("Threads (--jobs)", controls.jobs);

    controls.verbose = new QSpinBox(this);
    controls.verbose->setRange(0, 4);
    controls.verbose->setValue(3);
    layout->addRow("Verbose (--verbose)", controls.verbose);

    if (include_database)
    {
        layout->addRow("Database (--database)", BuildPathSelector(
            controls.database_path, controls.database_browse));
        controls.database_path->setText(
            QString::fromStdString(rhbm_gem::GetDefaultDatabasePath().string()));
    }

    layout->addRow("Output Folder (--folder)", BuildPathSelector(
        controls.folder_path, controls.folder_browse));
}

void MainWindow::ConnectUi()
{
    connect(m_command_list, &QListWidget::currentRowChanged, this, [this](int index)
    {
        m_command_stack->setCurrentIndex(index);
    });

    connect(m_run_button, &QPushButton::clicked, this, [this]()
    {
        StartExecution();
    });

    m_execution_poll_timer = new QTimer(this);
    m_execution_poll_timer->setInterval(50);
    connect(m_execution_poll_timer, &QTimer::timeout, this, [this]()
    {
        PollExecutionState();
    });

    connect(
        m_map_simulation.model_browse,
        &QPushButton::clicked,
        this,
        [this]() { BrowseForInputFile(m_map_simulation.model_path, "Select model file"); });
    connect(
        m_map_simulation.common.folder_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForDirectory(
                m_map_simulation.common.folder_path,
                "Select map simulation output folder");
        });

    connect(
        m_potential_analysis.model_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForInputFile(m_potential_analysis.model_path, "Select model file");
        });
    connect(
        m_potential_analysis.map_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForInputFile(m_potential_analysis.map_path, "Select map file");
        });
    connect(
        m_potential_analysis.training_report_dir_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForDirectory(
                m_potential_analysis.training_report_dir,
                "Select alpha training report directory");
        });
    connect(
        m_potential_analysis.common.database_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForOutputFile(
                m_potential_analysis.common.database_path,
                "Select output database path");
        });
    connect(
        m_potential_analysis.common.folder_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForDirectory(
                m_potential_analysis.common.folder_path,
                "Select potential analysis output folder");
        });
    connect(
        m_potential_analysis.simulation_flag,
        &QCheckBox::toggled,
        this,
        [this](bool enabled)
        {
            m_potential_analysis.simulated_map_resolution->setEnabled(enabled);
        });

    connect(
        m_result_dump.common.database_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForOutputFile(
                m_result_dump.common.database_path,
                "Select output database path");
        });
    connect(
        m_result_dump.common.folder_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForDirectory(
                m_result_dump.common.folder_path,
                "Select result dump output folder");
        });
    connect(
        m_result_dump.map_browse,
        &QPushButton::clicked,
        this,
        [this]()
        {
            BrowseForInputFile(m_result_dump.map_path, "Select map file");
        });
}

void MainWindow::BrowseForInputFile(QLineEdit * target, const QString & caption)
{
    const QString selected_path{
        QFileDialog::getOpenFileName(this, caption, target->text().trimmed())
    };
    if (!selected_path.isEmpty())
    {
        target->setText(selected_path);
    }
}

void MainWindow::BrowseForOutputFile(QLineEdit * target, const QString & caption)
{
    const QString selected_path{
        QFileDialog::getSaveFileName(this, caption, target->text().trimmed())
    };
    if (!selected_path.isEmpty())
    {
        target->setText(selected_path);
    }
}

void MainWindow::BrowseForDirectory(QLineEdit * target, const QString & caption)
{
    const QString selected_path{
        QFileDialog::getExistingDirectory(this, caption, target->text().trimmed())
    };
    if (!selected_path.isEmpty())
    {
        target->setText(selected_path);
    }
}

void MainWindow::StartExecution()
{
    if (m_execution_running) return;

    const int command_index{ m_command_stack->currentIndex() };
    if (command_index < 0 || static_cast<std::size_t>(command_index) >= m_gui_commands.size())
    {
        return;
    }
    const auto & command{ m_gui_commands[static_cast<std::size_t>(command_index)] };
    m_active_command_name = CommandName(command.id);
    m_execution_future = std::async(std::launch::async, [runner = command.run]()
    {
        return runner();
    });

    m_execution_running = true;
    m_execution_poll_timer->start();
    m_run_button->setEnabled(false);
    m_command_list->setEnabled(false);
    m_status_label->setText("Running " + m_active_command_name + " ...");
    const QString timestamp{ QDateTime::currentDateTime().toString(Qt::ISODate) };
    m_result_view->appendPlainText("[" + timestamp + "] Start: " + m_active_command_name);
}

void MainWindow::PollExecutionState()
{
    if (!m_execution_running || !m_execution_future.valid())
    {
        return;
    }
    const auto status{ m_execution_future.wait_for(std::chrono::milliseconds{ 0 }) };
    if (status != std::future_status::ready)
    {
        return;
    }

    m_execution_poll_timer->stop();
    m_execution_running = false;
    OnExecutionFinished(m_execution_future.get());
}

void MainWindow::OnExecutionFinished(const rhbm_gem::ExecutionReport & result)
{
    const QString summary{ FormatExecutionSummary(result) };
    const QString timestamp{ QDateTime::currentDateTime().toString(Qt::ISODate) };
    m_result_view->appendPlainText(
        "[" + timestamp + "] Finished: " + m_active_command_name + "\n" + summary + "\n");

    if (result.prepared && result.executed)
    {
        m_status_label->setText("Completed: " + m_active_command_name);
    }
    else if (!result.prepared)
    {
        m_status_label->setText("Validation failed: " + m_active_command_name);
    }
    else
    {
        m_status_label->setText("Execution failed: " + m_active_command_name);
    }

    m_command_list->setEnabled(true);
    m_run_button->setEnabled(true);
}

rhbm_gem::MapSimulationRequest MainWindow::BuildMapSimulationRequest() const
{
    rhbm_gem::MapSimulationRequest request;
    request.common.thread_size = m_map_simulation.common.jobs->value();
    request.common.verbose_level = m_map_simulation.common.verbose->value();
    request.common.folder_path = ReadPath(m_map_simulation.common.folder_path);
    request.model_file_path = ReadPath(m_map_simulation.model_path);
    request.map_file_name = ReadString(m_map_simulation.map_name);
    request.potential_model_choice = CurrentEnumValue<rhbm_gem::PotentialModel>(
        m_map_simulation.potential_model);
    request.partial_charge_choice = CurrentEnumValue<rhbm_gem::PartialCharge>(
        m_map_simulation.partial_charge);
    request.cutoff_distance = m_map_simulation.cutoff_distance->value();
    request.grid_spacing = m_map_simulation.grid_spacing->value();
    request.blurring_width_list = ReadString(m_map_simulation.blurring_widths);
    return request;
}

rhbm_gem::PotentialAnalysisRequest MainWindow::BuildPotentialAnalysisRequest() const
{
    rhbm_gem::PotentialAnalysisRequest request;
    request.common.thread_size = m_potential_analysis.common.jobs->value();
    request.common.verbose_level = m_potential_analysis.common.verbose->value();
    request.common.database_path = ReadPath(m_potential_analysis.common.database_path);
    request.common.folder_path = ReadPath(m_potential_analysis.common.folder_path);

    request.model_file_path = ReadPath(m_potential_analysis.model_path);
    request.map_file_path = ReadPath(m_potential_analysis.map_path);
    request.saved_key_tag = ReadString(m_potential_analysis.saved_key_tag);
    request.training_report_dir = ReadPath(m_potential_analysis.training_report_dir);
    request.simulation_flag = m_potential_analysis.simulation_flag->isChecked();
    request.simulated_map_resolution = m_potential_analysis.simulated_map_resolution->value();
    request.training_alpha_flag = m_potential_analysis.training_alpha_flag->isChecked();
    request.asymmetry_flag = m_potential_analysis.asymmetry_flag->isChecked();
    request.sampling_size = m_potential_analysis.sampling_size->value();
    request.sampling_range_min = m_potential_analysis.sampling_range_min->value();
    request.sampling_range_max = m_potential_analysis.sampling_range_max->value();
    request.sampling_height = m_potential_analysis.sampling_height->value();
    request.fit_range_min = m_potential_analysis.fit_range_min->value();
    request.fit_range_max = m_potential_analysis.fit_range_max->value();
    request.alpha_r = m_potential_analysis.alpha_r->value();
    request.alpha_g = m_potential_analysis.alpha_g->value();
    return request;
}

rhbm_gem::ResultDumpRequest MainWindow::BuildResultDumpRequest() const
{
    rhbm_gem::ResultDumpRequest request;
    request.common.thread_size = m_result_dump.common.jobs->value();
    request.common.verbose_level = m_result_dump.common.verbose->value();
    request.common.database_path = ReadPath(m_result_dump.common.database_path);
    request.common.folder_path = ReadPath(m_result_dump.common.folder_path);

    request.printer_choice = CurrentEnumValue<rhbm_gem::PrinterType>(m_result_dump.printer_choice);
    request.model_key_tag_list = ReadString(m_result_dump.model_key_list);
    request.map_file_path = ReadPath(m_result_dump.map_path);
    return request;
}

QString MainWindow::FormatExecutionSummary(const rhbm_gem::ExecutionReport & result)
{
    QString summary;
    summary += QString("Prepared: %1\n").arg(result.prepared ? "true" : "false");
    summary += QString("Executed: %1\n").arg(result.executed ? "true" : "false");
    summary += "Validation Issues:\n";
    summary += FormatValidationIssues(result.validation_issues);
    return summary;
}

QString MainWindow::FormatValidationIssues(
    const std::vector<rhbm_gem::ValidationIssue> & issues)
{
    if (issues.empty())
    {
        return "  (none)";
    }

    QString lines;
    for (const auto & issue : issues)
    {
        lines += "  - [";
        lines += ValidationPhaseLabel(issue.phase);
        lines += "/";
        lines += ValidationLevelLabel(issue.level);
        lines += "] ";
        lines += QString::fromStdString(issue.option_name);
        lines += ": ";
        lines += QString::fromStdString(issue.message);
        if (issue.auto_corrected)
        {
            lines += " (auto-corrected)";
        }
        lines += "\n";
    }

    lines.chop(1);
    return lines;
}

QString MainWindow::ValidationPhaseLabel(rhbm_gem::ValidationPhase phase)
{
    switch (phase)
    {
    case rhbm_gem::ValidationPhase::Parse:
        return "parse";
    case rhbm_gem::ValidationPhase::Prepare:
        return "prepare";
    default:
        return "unknown";
    }
}

QString MainWindow::ValidationLevelLabel(LogLevel level)
{
    switch (level)
    {
    case LogLevel::Error:
        return "error";
    case LogLevel::Warning:
        return "warning";
    case LogLevel::Notice:
        return "notice";
    case LogLevel::Info:
        return "info";
    case LogLevel::Debug:
        return "debug";
    default:
        return "unknown";
    }
}
