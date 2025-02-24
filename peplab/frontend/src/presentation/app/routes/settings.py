from flask import request, session, jsonify, render_template, current_app, flash
from flask import Blueprint
from flask_wtf.csrf import CSRFProtect

bp = Blueprint("settings", __name__)


@bp.route("/save_appearance", methods=["POST"])
def save_appearance():
    try:
        theme = request.form.get("theme", "system")
        if theme not in ["system", "light", "dark"]:
            return (
                jsonify({"status": "error", "message": "Invalid theme selection"}),
                400,
            )

        session["theme"] = theme
        return jsonify(
            {"status": "success", "theme": theme, "message": "Theme saved successfully"}
        )
    except Exception as e:
        current_app.logger.error(f"Error saving theme: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500


@bp.route("/settings")
def settings():
    current_theme = session.get("theme", "system")
    return render_template("settings/settings.html", current_theme=current_theme)
