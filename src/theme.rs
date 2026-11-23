use ratatui::style::{Color, Modifier, Style};
use ratatui::text::Span;
use ratatui::widgets::{Block, BorderType, Borders};

pub const BG_TRANSPARENT: Color = Color::Reset;
pub const SURFACE_DARK: Color = Color::Rgb(35, 42, 46);
pub const SURFACE_ALT: Color = Color::Rgb(45, 53, 59);
pub const BORDER_LEFT: Color = Color::Rgb(66, 88, 66);
pub const BORDER_RIGHT: Color = Color::Rgb(232, 176, 96);

pub const ACCENT_GREEN: Color = Color::Rgb(132, 204, 122);
pub const ACCENT_AMBER: Color = Color::Rgb(232, 176, 96);
pub const ACCENT_AQUA: Color = Color::Rgb(96, 176, 144);
pub const ACCENT_YELLOW: Color = Color::Rgb(224, 200, 120);
pub const ACCENT_CORAL: Color = Color::Rgb(208, 112, 112);
pub const ACCENT_PURPLE: Color = Color::Rgb(176, 128, 176);

pub const TEXT_PRIMARY: Color = Color::Rgb(234, 238, 223);
pub const TEXT_MUTED: Color = Color::Rgb(120, 136, 120);
pub const DARK_TEXT: Color = Color::Rgb(20, 28, 22);

pub const TRACE_COLORS: [Color; 6] = [
    ACCENT_GREEN,
    ACCENT_AMBER,
    ACCENT_AQUA,
    ACCENT_PURPLE,
    ACCENT_CORAL,
    ACCENT_YELLOW,
];

pub fn trace_color(idx: usize) -> Color {
    TRACE_COLORS[idx % TRACE_COLORS.len()]
}

pub fn panel_block<'a>(title: impl Into<String>, border_color: Color) -> Block<'a> {
    Block::default()
        .borders(Borders::ALL)
        .border_type(BorderType::Plain)
        .border_style(Style::default().fg(border_color))
        .title(Span::styled(
            format!(" {} ", title.into()),
            Style::default()
                .fg(border_color)
                .add_modifier(Modifier::BOLD),
        ))
}

pub fn modal_block<'a>(title: impl Into<String>, border_color: Color) -> Block<'a> {
    Block::default()
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .style(Style::default().bg(SURFACE_DARK).fg(TEXT_PRIMARY))
        .border_style(Style::default().fg(border_color))
        .title(Span::styled(
            format!(" {} ", title.into()),
            Style::default()
                .fg(border_color)
                .bg(SURFACE_DARK)
                .add_modifier(Modifier::BOLD),
        ))
}

pub fn tab_span<'a>(label: &'a str, active: bool) -> Span<'a> {
    if active {
        Span::styled(
            format!(" {} ", label),
            Style::default()
                .fg(DARK_TEXT)
                .bg(ACCENT_GREEN)
                .add_modifier(Modifier::BOLD),
        )
    } else {
        Span::styled(
            format!(" {} ", label),
            Style::default().fg(TEXT_MUTED).bg(SURFACE_DARK),
        )
    }
}

pub fn key_badge<'a>(label: &'a str, fg: Color) -> Span<'a> {
    Span::styled(
        format!(" {} ", label),
        Style::default()
            .fg(fg)
            .bg(SURFACE_DARK)
            .add_modifier(Modifier::BOLD),
    )
}
